from graphviz import Digraph
import itertools
from math import log2
from pysat.solvers import Glucose3
from pysat.card import *

class SynthesisException(Exception):

    def __init__(self, message):
        self.message = message

class node:
    '''
    A generic node class, used by both clocking scheme graphs
    and logic networks.
    '''
    def __init__(self, *, coords=None, is_pi=False, is_po=False):
        self.coords = coords
        self.is_pi = is_pi
        self.virtual_fanin = []
        self.virtual_fanout = []
        self.is_border_node = False
        self.is_po = is_po
        self.fanin = {}
        self.is_designated_pi = False
        self.is_designated_po = False
        self.gate_type = None

    def __repr__(self):
        if self.is_pi:
            return 'PI{}'.format(self.coords)
        else:
            return '<node {}>'.format(self.coords)
        
    def __lt__(self, other): 
        return self.coords < other.coords

class logic_network:
    '''
    A logic_network is the result of synthesis. Its design is similar
    to a clocking scheme graph. A major difference is that it cannot
    contain cycles. However, it can also be accessed using the same
    tile coordinate based API.
    '''
    
    def __init__(self, shape, nr_pis, nr_pos):
        '''
        Creates a new logic network according to the following
        specifications:

        shape: a 2-tuple containing the size of the grid
               (width x height)
        nr_pis: number of PIs
        nr_pos: number of POs

        Furthermore, it has the following properties:
        self.nodes: A list of all the nodes in the logic 
                    network, including the PIs.
        self.node_map: A map from tile coordinates to nodes
                       nodes in the logic network. E.g. to
                       access the node corresponding to
                       tile (0,0) one refers to 
                       node_map[(0,0)].
        self.po_map: A list of size nr_pos mapping output
                     indices to nodes in the network. E.g.
                     to access the first output, one refers
                     to po_map[0].
        '''
        self.shape = shape
        self.nr_pis = nr_pis
        self.nr_pos = nr_pos
        
        self.nodes = [node(coords=i, is_pi=True) for i in range(nr_pis)] 
        self.node_map = {}
        for y in range(shape[1]):
            for x in range(shape[0]):
                n = node(coords=(x, y))
                if x == 0 or x == (shape[0] - 1) or y == 0 or y == (shape[1] - 1):
                    n.is_border_node = True
                self.nodes.append(n)
                self.node_map[(x,y)] = n

        self.po_map = [None] * nr_pos

    def set_output(self, h, coords):
        '''
        Marks the node at coords as the h-th output of the network.
        '''
        assert(h < self.nr_pos)
        n = self.node_map[coords]
        n.is_po = True
        self.po_map[h] = n

    def __repr__(self):
        r = '\n'
        for n in self.nodes:
            if n.is_pi:
                continue
            r += '<node {}, function={}'.format(n.coords, n.function)
            if 0 in n.fanin:
                r += ', fanin[0]={}'.format(n.fanin[0])
            if 1 in n.fanin:
                r += ', fanin[1]={}'.format(n.fanin[1])
            if 2 in n.fanin:
                r += ', fanin[2]={}'.format(n.fanin[2])
            r += '>\n'
        return r

    def has_border_io(self):
        '''
        Checks if only border nodes are connected to PIs and POs. Returns
        True of this is the case and False otherwise.
        '''
        for n in self.nodes:
            if n.is_pi:
                continue
            if not n.is_border_node:
                for innode in n.fanin.values():
                    if innode.is_pi:
                        return False
        for n in self.po_map:
            if not n.is_border_node:
                return False
        return True

    def to_png(self, filename):
        '''
        Creates a PNG of the logic network using Graphviz.
        In the resulting PNG, all border nodes are filled in
        with a gray color. All internal nodes are white.
        PO nodes are marked by a double border. Every non-PI
        node is also marked with its tile-space coordinates
        as well as the function it computes.
        '''
        dot = Digraph()
        dot.attr('node', shape='circle')
        
        pi_counter = 0
        for n in self.nodes:
            if n.is_pi:
                continue
            name = 'N_{}_{}'.format(n.coords[0], n.coords[1])
            label = '{}\n{}\n{}'.format(n.coords, n.gate_type, str(n.function))
            fill = 'gray' if n.is_border_node else 'white'
            if n.is_po:
                dot.node(name, label, shape='doublecircle', fillcolor=fill, style='filled')
            else:
                dot.node(name, label, fillcolor=fill, style='filled')
            for innode in n.fanin.values():
                if innode.is_pi:
                    inname = 'PI_{}'.format(pi_counter)
                    pi_counter += 1
                    inlabel = 'PI{}'.format(innode.coords)
                    dot.node(inname, inlabel)
                    dot.edge(inname, name)
                else:
                    inname = 'N_{}_{}'.format(innode.coords[0], innode.coords[1])
                    dot.edge(inname, name)

        dot.render(filename=filename, format='png', cleanup=True)

class scheme_graph:
    '''
    A scheme_graph (short for clocking-scheme graph) is used
    to specify a clocking scheme and to synthesize logic
    networks according to that specification.
    '''

    # Mapping gate types to their local truth tables
    gate_tts = {
        'NOT':  [1, 0, 1, 0, 1, 0, 1, 0],
        'AND':  [0, 0, 0, 1, 0, 0, 0, 1],
        'OR':   [0, 1, 1, 1, 0, 1, 1, 1],
        'MAJ':  [0, 0, 0, 1, 0, 1, 1, 1],
        'WIRE': [0, 1, 0, 1, 0, 1, 0, 1]
    }

    # Reverse map of truth tables to gate types.
    tt_gate_type = {
        (1, 0, 1, 0, 1, 0, 1, 0): 'NOT',
        (0, 0, 0, 1, 0, 0, 0, 1): 'AND',
        (0, 1, 1, 1, 0, 1, 1, 1): 'OR',
        (0, 0, 0, 1, 0, 1, 1, 1): 'MAJ',
        (0, 1, 0, 1, 0, 1, 0, 1): 'WIRE'
    }

    # Mapping gate types to the corresponding number of fanins.
    gate_fanin_range = {
        'NOT': 1,
        'AND': 2,
        'OR': 2,
        'MAJ': 3,
        'WIRE': 1
    }

    def __init__(self, *, shape=(1,1), nr_pis=1, nr_pos=1, border_io=False,
                 enable_not=True, enable_and=True, enable_or=True,
                 enable_maj=True, enable_crossings=False):
        '''
        Creates a new clocking scheme graph according to specifications.
        Defines the following properties.

        shape: A 2-tuple specifying the dimensions of the clocking scheme.
        nr_pis: Number of PIs.
        nr_pos: Number of POs.
        border_io: True iff only border nodes can have PI fanin.
        enable_not: Enable synthesis of NOT gates.
        enable_and: Enable synthesis of AND gates.
        enable_or: Enable synthesis of OR gates.
        enable_maj: Enable synthesis of MAJ gates.
        enable_crossings: Enable wire crossings.
        '''
        self.shape = shape
        self.nr_pis = nr_pis
        self.nr_pos = nr_pos
        
        self.nodes = [node(coords=i,is_pi=True) for i in range(nr_pis)] 
        self.node_map = {}
        for y in range(shape[1]):
            for x in range(shape[0]):
                n = node(coords=(x, y))
                if x == 0 or x == (shape[0] - 1) or y == 0 or y == (shape[1] - 1):
                    n.is_border_node = True
                self.nodes.append(n)
                self.node_map[(x,y)] = n

        self.border_io = border_io
        self.enable_not = enable_not
        self.enable_and = enable_and
        self.enable_or = enable_or
        self.enable_maj = enable_maj
        self.enable_crossings = enable_crossings
    
    def add_virtual_edge(self, coords1, coords2):
        '''
        Adds a virtual edge from the node corresponding to the tile at
        coords1 to the node corresponding to the tile at coords2.
        A virtual edge specifies that the node at coords2 may have
        the node at coords1 in its fanin. However, it does not force
        this to happen. Hence, the connection is virtual and may be
        actualized by the synthesis process.
        '''
        node1 = self.node_map[coords1]
        node2 = self.node_map[coords2]
        node1.virtual_fanout.append(node2)
        node2.virtual_fanin.append(node1)

    def designate_pi(self, coords):
        '''
        Designates the node at the given tile coordinates as an PI
        pin. It will be hardcoded to act as a wire element which is
        connected to a PI.
        '''
        n = self.node_map[coords]
        n.is_designated_pi = True

    def designate_po(self, coords):
        '''
        Designates the node at the given tile coordinates as a PO
        pint. It will be hardcoded to act as a wire element which
        connected to a PO.
        '''
        n = self.node_map[coords]
        n.is_designated_po = True

    def to_png(self, filename):
        '''
        Creates a PNG of the graph underlying the clock scheme 
        using Graphviz.
        '''
        dot = Digraph()
        dot.attr('node', shape='box')
        dot.attr(splines='ortho')

        for y in range(self.shape[1]):
            with dot.subgraph() as s:
                s.attr(rank='same')
                for x in range(self.shape[0]):
                    n = self.node_map[(x, y)]
                    name = 'N_{}_{}'.format(x, y)
                    label = str((x, y))
                    s.node(name, label)
                    if x > 0:
                        prevname = 'N_{}_{}'.format(x-1, y)
                        s.edge(prevname, name, style='invis')


        for y in range(self.shape[1]):
            for x in range(self.shape[0]):
                n1 = self.node_map[(x, y)]
                n1name = 'N_{}_{}'.format(x, y)
                for n2 in n1.virtual_fanout:
                    n2name = 'N_{}_{}'.format(n2.coords[0], n2.coords[1])
                    dot.edge(n1name, n2name)

        dot.render(filename=filename, format='png', cleanup=True)

    def synthesize(self, functions, verbosity=0):
        '''
        Synthesizes a logic network according to the clocking scheme
        specifications encoded in the graph and the functional
        specification encoded by the truth tables in the functions
        list.
        
        NOTE: this function may be called multiple times, which
        will result in it generating zero or more logic networks.
        '''
        assert(len(functions) > 0)
        assert(len(functions[0]) % 2 == 0)
        assert(round(log2(len(functions[0]))) == self.nr_pis)

        sim_vars = {}
        input_sim_vars = {}
        lut_vars = {}
        var_idx = 1

        enabled_gates = ['WIRE']
        if self.enable_not:
            enabled_gates.append('NOT')
        if self.enable_and:
            enabled_gates.append('AND')
        if self.enable_or:
            enabled_gates.append('OR')
        if self.enable_maj:
            enabled_gates.append('MAJ')
        assert(len(enabled_gates) > 1)
        nr_gate_types = len(enabled_gates)
        if verbosity > 0:
            print('enabled_gates={}'.format(enabled_gates))

        # Pre-process specified functions, make sure truth tables are
        # at least size 8.
        if len(functions[0]) == 4:
            for i in range(len(functions)):
                functions[i] = functions[i] + functions[i]

        # Create simulation, input simulation, and LUT variables
        nr_sim_vars = len(functions[0])
        nr_input_sim_vars = len(functions[0]) * 3
        nr_lut_vars = 8
        for n in self.nodes:
            nsimvars = [0] * nr_sim_vars
            for i in range(nr_sim_vars):
                nsimvars[i] = var_idx
                var_idx += 1
            sim_vars[n] = nsimvars

            if n.is_pi:
                continue

            ninputsimvars = {}
            for i in range(3):
                nsimvars = [0] * nr_sim_vars
                for j in range(nr_sim_vars):
                    nsimvars[j] = var_idx
                    var_idx += 1
                ninputsimvars[i] = nsimvars
            input_sim_vars[n] = ninputsimvars

            nlutvars = [0] * nr_lut_vars
            for i in range(nr_lut_vars):
                nlutvars[i] = var_idx
                var_idx += 1
            lut_vars[n] = nlutvars

        nr_sim_and_lut_vars = var_idx - 1
        if verbosity > 0:
            print('Created {} SIM and LUT vars'.format(nr_sim_and_lut_vars))

        # Create selection variables
        svar_map = {} # Maps variables to PIs or internal node_map
        fanin_options = {}
        pi_fanin_options = [self.nodes[x] for x in range(self.nr_pis)]
        for n in self.nodes:
            if n.is_pi:
                continue
            fanin_options[n] = {}
            svar_map[n] = {}
            for k in range(3):
                if n.is_designated_pi:
                    fanin_options[n][k] = pi_fanin_options
                elif (self.border_io and n.is_border_node) or not self.border_io:
                    fanin_options[n][k] = pi_fanin_options + n.virtual_fanin
                else:
                    fanin_options[n][k] = n.virtual_fanin
                nr_fanin_options = len(fanin_options[n][k])
                if verbosity > 0:
                    print('fanin options for {}[{}]: {}'.format(n.coords, k, fanin_options[n][k]))
                    print('nr fanin options: {}'.format(nr_fanin_options))
                assert(nr_fanin_options > 0)
                # Create variables for the 3 potentially used fanin ports
                svar_map[n][k] = [0] * nr_fanin_options
                for i in range(nr_fanin_options):
                    svar_map[n][k][i] = var_idx
                    var_idx += 1

        nr_selection_vars = var_idx - 1 - nr_sim_and_lut_vars
        if verbosity > 0:
            print('Created {} selection vars'.format(nr_selection_vars))

        # Create output variables
        nr_outputs = len(functions)
        out_vars = {}
        for h in range(nr_outputs):
            houtvars = {}
            for y in range(self.shape[1]):
                for x in range(self.shape[0]):
                    n = self.node_map[(x,y)]
                    if self.border_io and not n.is_border_node:
                        continue
                    houtvars[(x,y)] = var_idx
                    var_idx += 1
            out_vars[h] = houtvars

        nr_out_vars = var_idx - 1 - nr_sim_and_lut_vars - nr_selection_vars
        if verbosity > 0:
            print('Created {} output vars'.format(nr_out_vars))

        gate_type_vars = {}
        for n in self.nodes:
            if n.is_pi:
                continue
            ngatetypevars = [0] * nr_gate_types
            for i in range(nr_gate_types):
                ngatetypevars[i] = var_idx
                var_idx += 1
            gate_type_vars[n] = ngatetypevars

        nr_gate_type_vars = var_idx - 1 - nr_sim_and_lut_vars - nr_selection_vars - nr_out_vars
        if verbosity > 0:
            print('Created {} gate-type vars'.format(nr_gate_type_vars))

        # Create graph connection and path variables
        path_vars = {}
        for n in self.nodes:
            if n.is_pi:
                continue
            path_vars[n] = {}
            for np in self.nodes:
                if np.is_pi:
                    continue
                path_vars[n][np] = var_idx
                var_idx += 1

        nr_path_vars = var_idx - 1 - nr_sim_and_lut_vars - nr_selection_vars - nr_out_vars - nr_gate_type_vars
        if verbosity > 0:
            print('Created {} path vars'.format(nr_path_vars))

        clauses = []
        # Create simulation propagation constraints
        for n in self.nodes:
            if n.is_pi:
                continue
            for tt_idx in range(nr_sim_vars):
                permutations = list(itertools.product('01', repeat=4))
                for permutation in permutations:
                    function_output = int(permutation[0])
                    const_vals = []
                    f_idx = 0
                    for i in range(3):
                        const_val = int(permutation[i+1])
                        const_vals.append(const_val)
                        f_idx += (const_val << i)

                    clause = []
                    if function_output == 1:
                        clause.append(-sim_vars[n][tt_idx])
                    else:
                        clause.append(sim_vars[n][tt_idx])
                    for i in range(len(const_vals)):
                        if const_vals[i] == 1:
                            clause.append(-input_sim_vars[n][i][tt_idx])
                        else:
                            clause.append(input_sim_vars[n][i][tt_idx])
                    if function_output == 1:
                        clause.append(lut_vars[n][f_idx])
                    else:
                        clause.append(-lut_vars[n][f_idx])
                    clauses.append(clause)

        # Make sure each node selects some fanin on each of its input ports
        for n in self.nodes:
            if n.is_pi:
                continue
            for k in range(3):
                clauses.append(svar_map[n][k])

        # Add fanin selection constraints
        for y in range(self.shape[1]):
            for x in range(self.shape[0]):
                n = self.node_map[(x,y)]
                for k in range(3):
                    options = fanin_options[n][k]
                    nr_fanin_options = len(options)
                    for i in range(nr_fanin_options):
                        innode = options[i]
                        # Suppose fanin k of node n selects option i
                        # Then, for all tt_idx, input truth table of node n
                        # at tt_idx must equal truth table of i at tt_idx
                        for tt_idx in range(nr_sim_vars):
                            clause = [0] * 3
                            clause[0] = -svar_map[n][k][i]
                            clause[1] = -input_sim_vars[n][k][tt_idx]
                            clause[2] = sim_vars[innode][tt_idx]
                            clauses.append(clause)
                            clause = [0] * 3
                            clause[0] = -svar_map[n][k][i]
                            clause[1] = input_sim_vars[n][k][tt_idx]
                            clause[2] = -sim_vars[innode][tt_idx]
                            clauses.append(clause)

        # Create cycle-prevention constraints
        for n in self.nodes:
            if n.is_pi:
                continue
            for k in range(3):
                options = fanin_options[n][k]
                nr_fanin_options = len(options)
                for i in range(nr_fanin_options):
                    innode = options[i]
                    if innode.is_pi:
                        continue
                    # Suppose fanin k of node n selects option i.
                    # Then, there is a connection from node i to node
                    # n.  Therefore, we must set the variable encoding
                    # path[i][n] to true.
                    clause = [0] * 2
                    clause[0] = -svar_map[n][k][i]
                    clause[1] = path_vars[innode][n]
                    clauses.append(clause)

        # For all nodes n, if is path[n][n'] is true and
        # path[n'][n''] is true, then path[n][n''] is true as
        # well.
        for n in self.nodes:
            if n.is_pi:
                continue
            for np in self.nodes:
                if np.is_pi:
                    continue
                for npp in self.nodes:
                    if npp.is_pi:
                        continue
                    clauses.append([-path_vars[n][np], -path_vars[np][npp], path_vars[n][npp]])
            # Ensure there is no path from n to itself.
            clauses.append([-path_vars[n][n]])

        # Fix input vars
        for var in range(self.nr_pis):
            for idx in range(nr_sim_vars):
                if idx & (1 << var):
                    clauses.append([sim_vars[self.nodes[var]][idx]])
                else:
                    clauses.append([-sim_vars[self.nodes[var]][idx]])

        # Fix output vars
        for h in range(nr_outputs):
            houtvars = out_vars[h]
            # Ensure that output h points to exactly one gate
            cnf = CardEnc.equals(lits=list(houtvars.values()), encoding=EncType.pairwise)
            for clause in cnf.clauses:
                clauses.append(clause)
            # Encode the different choices
            for y in range(self.shape[1]):
                for x in range(self.shape[0]):
                    n = self.node_map[(x,y)]
                    if self.border_io and not n.is_border_node:
                        continue
                    # If output h points to gate (x,y), then the truth table
                    # of (x,y) must agree with that of function h
                    houtvar = houtvars[(x,y)]
                    for idx in range(nr_sim_vars):
                        if functions[h][idx] == 1:
                            clauses.append([-houtvar, sim_vars[n][idx]])
                        else:
                            clauses.append([-houtvar, -sim_vars[n][idx]])    

        # Add gate constraints
        for n in self.nodes:
            if n.is_pi:
                continue
            # Gate n must pick exactly one of the enabled gate types
            ngatetypevars = gate_type_vars[n]
            cnf = CardEnc.equals(lits=ngatetypevars, encoding=EncType.pairwise)
            for clause in cnf.clauses:
                clauses.append(clause)
            for i in range(nr_gate_types):
                gate_type = enabled_gates[i]
                gate_tt = self.gate_tts[gate_type]
                nlutvars = lut_vars[n]
                for j in range(8):
                    if gate_tt[j] == 1:
                        clauses.append([-ngatetypevars[i], nlutvars[j]])
                    else:
                        clauses.append([-ngatetypevars[i], -nlutvars[j]])

        # Add cardinality constraints on gate fanouts. These
        # constraints depend on the gate type. AND/OR/NOT/MAJ gates
        # are restricted to single fanout, while wires may have fanout
        # up to three.
        for n in self.nodes:
            if n.is_pi:
                continue
            fanout_vars = []
            for np in self.nodes:
                if np.is_pi or np == n:
                    continue
                for k in range(3):
                    options = fanin_options[np][k]
                    nr_fanin_options = len(options)
                    for i in range(nr_fanin_options):
                        if options[i] == n:
                            fanout_vars.append(svar_map[np][k][i])
            if verbosity > 1:
                print('fanout_vars {}: {}'.format(n.coords, fanout_vars))
            # Create cardinality constraints based on gate type.  The
            # first gate type var is always a WIRE.
            ngatetypevars = gate_type_vars[n]
            wire_var = ngatetypevars[0]
            cnf = CardEnc.atmost(lits=fanout_vars, encoding=EncType.pairwise, bound=3)
            for clause in cnf.clauses:
                clauses.append([-wire_var] + clause)
            for gt_var in ngatetypevars[1:]:
                cnf = CardEnc.atmost(lits=fanout_vars, encoding=EncType.pairwise)
                for clause in cnf.clauses:
                    clauses.append([-gt_var] + clause)

        # Create clauses to fix designated I/O pins
        for n in self.nodes:
            if n.is_designated_pi or n.is_designated_po:
                # The first gate type var always refers to WIRE and n
                # must be a WIRE element.
                ngatetypevars = gate_type_vars[n]
                clauses.append([ngatetypevars[0]])
            if n.is_designated_po:
                # At least one of the outputs must point to n.
                outpointers = []
                for h in range(nr_outputs):
                    houtvars = out_vars[h]
                    houtvar = houtvars[(x,y)]
                    outpointers.append(houtvar)
                clauses.append(outpointers)
                
        nr_clauses = len(clauses)
        if verbosity > 0:
            print('nr clauses: {}'.format(nr_clauses))

        solver = Glucose3()
        for clause in clauses:
            solver.add_clause(clause)

#        print(clauses)
        for model in solver.enum_models():
            # Decode network from solution
#            print(model)
            if verbosity > 1:
                for i in range(len(self.nodes)):
                    n = self.nodes[i]
                    if n.is_pi:
                        continue
                    print('tt {}: '.format(n.coords), end='')
                    for j in range(nr_sim_vars):
                        val = 1 if model[sim_vars[n][j]-1] > 0 else 0
                        print('{}'.format(val), end='')
                    print()
            
            net = logic_network(self.shape, self.nr_pis, len(functions))
            for h in range(nr_outputs):
                houtvars = out_vars[h]
                out_found = 0
                for y in range(self.shape[1]):
                    for x in range(self.shape[0]):
                        n = self.node_map[(x,y)]
                        if self.border_io and not n.is_border_node:
                            continue
                        houtvar = houtvars[(x,y)]
                        if model[houtvar-1] > 0:
                            if verbosity > 1:
                                print('out[{}] -> ({},{})'.format(h, x, y))
                            out_found += 1
                            net.set_output(h, (x,y))
                            
                # Every output must point somewhere
                assert(out_found == 1)

            for n in self.nodes:
                if n.is_pi:
                    continue
                lut_tt = [0] * 8
                for i in range(nr_lut_vars):
                    var = lut_vars[n][i]
                    model_val = model[var - 1]
                    if model_val > 0:
                        lut_tt[i] = 1
                    else:
                        lut_tt[i] = 0
                assert(tuple(lut_tt) in self.tt_gate_type)
                gate_type = self.tt_gate_type[tuple(lut_tt)]
                fanin_range = self.gate_fanin_range[gate_type]
                netnode = net.node_map[n.coords]
                netnode.is_border_node = n.is_border_node
                netnode.function = lut_tt
                netnode.gate_type = gate_type
                for k in range(fanin_range):
                    options = fanin_options[n][k]
                    for j in range(len(options)):
                        var = svar_map[n][k][j]
                        model_val = model[var - 1]
                        if model_val > 0:
                            innode = options[j]
                            if innode.is_pi:
                                netnode.fanin[k] = net.nodes[innode.coords]
                            else:
                                netnode.fanin[k] = net.node_map[innode.coords]
                            if verbosity > 1:
                                print('{}[{}] = {} (var {})'.format(n.coords, k, innode.coords, var))
            yield net
