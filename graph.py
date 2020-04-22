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
        self.fanout = []
        self.gate_type = None

    def set_fanin(self, k, innode):
        '''
        Sets the k-th fanin of this node to innode and updates the fanout
        of innode by appending this node to it.
        '''
        self.fanin[k] = innode
        innode.fanout.append(self)

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
        True if this is the case and False otherwise.
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

    def has_designated_pi(self):
        '''
        Checks if only WIREs are connected to PIs. Moreover, verifies that
        those designated PI WIREs have only one fanout. Returns True
        of this is the case and False otherwise.
        '''
        for n in self.nodes:
            if n.is_pi:
                continue
            if n.gate_type != 'WIRE':
                for innode in n.fanin.values():
                    if innode.is_pi:
                        return False
            else:
                nr_fanout = len(n.fanout)
                for innode in n.fanin.values():
                    if innode.is_pi and nr_fanout > 1:
                        return False
        return True

    def has_designated_po(self):
        '''
        Checks if only WIREs are connected to POs. Moreover, verifies that
        those designated PO WIREs have no other fanout. Returns True
        of this is the case and False otherwise.
        '''
        for n in self.po_map:
            if n.gate_type != 'WIRE':
                return False
            if len(n.fanout) > 0:
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

    def __init__(self, *, shape=(1,1), border_io=False,
                 enable_wire=True, enable_not=True, enable_and=True,
                 enable_or=True, enable_maj=True,
                 enable_crossings=False, designated_pi=False,
                 designated_po=False):
        '''
        Creates a new clocking scheme graph according to specifications.
        Defines the following properties.

        shape: A 2-tuple specifying the dimensions of the clocking scheme.
        border_io: True iff only border nodes can have PI fanin.
        enable_not: Enable synthesis of WIREs.
        enable_not: Enable synthesis of NOT gates.
        enable_and: Enable synthesis of AND gates.
        enable_or: Enable synthesis of OR gates.
        enable_maj: Enable synthesis of MAJ gates.
        enable_crossings: Enable wire crossings.
        designated_pi: True iff only WIRES can have PI fanin.
        designated_po: True iff only WIRES can have PO fanout.
        '''
        self.shape = shape
        self.node_map = {}
        for y in range(shape[1]):
            for x in range(shape[0]):
                n = node(coords=(x, y))
                if x == 0 or x == (shape[0] - 1) or y == 0 or y == (shape[1] - 1):
                    n.is_border_node = True
                self.node_map[(x,y)] = n

        self.enable_wire = enable_wire
        self.border_io = border_io
        self.enable_not = enable_not
        self.enable_and = enable_and
        self.enable_or = enable_or
        self.enable_maj = enable_maj
        self.enable_crossings = enable_crossings
        self.designated_pi = designated_pi
        self.designated_po = designated_po
    
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

        self.nr_pis = round(log2(len(functions[0])))
        self.nr_pos = len(functions)
        
        self.nodes = [node(coords=i,is_pi=True) for i in range(self.nr_pis)] 
        for y in range(self.shape[1]):
            for x in range(self.shape[0]):
                self.nodes.append(self.node_map[(x,y)])

        sim_vars = {}
        lut_vars = {}
        rev_var_map = {}
        var_idx = 1

        enabled_gates = []
        if self.enable_wire:
            enabled_gates.append('WIRE')
        if self.enable_not:
            enabled_gates.append('NOT')
        if self.enable_and:
            enabled_gates.append('AND')
        if self.enable_or:
            enabled_gates.append('OR')
        if self.enable_maj:
            enabled_gates.append('MAJ')
        assert(len(enabled_gates) > 0)
        nr_gate_types = len(enabled_gates)
        # Determine the minimum and maximum number of fanins that are
        # enabled.
        gate_fanin_options = [self.gate_fanin_range[t] for t in enabled_gates]
        min_fanin = min(gate_fanin_options)
        max_fanin = max(gate_fanin_options)
        
        if verbosity > 0:
            print('enabled_gates={}'.format(enabled_gates))

        # Pre-process specified functions, make sure truth tables are
        # at least size 8.
        if len(functions[0]) == 4:
            for i in range(len(functions)):
                functions[i] = functions[i] + functions[i]

        # Create simulation, input simulation, and LUT variables
        nr_sim_vars = len(functions[0])
        nr_lut_vars = 8
        for n in self.nodes:
            nsimvars = [0] * nr_sim_vars
            for i in range(nr_sim_vars):
                nsimvars[i] = var_idx
                rev_var_map[var_idx] = 'simvars[{}][{}] = {}'.format(n.coords, i, var_idx)
                var_idx += 1
            sim_vars[n] = nsimvars

            if n.is_pi:
                continue

            nlutvars = [0] * nr_lut_vars
            for i in range(nr_lut_vars):
                nlutvars[i] = var_idx
                rev_var_map[var_idx] = 'nlutvars[{}] = {}'.format(i, var_idx)
                var_idx += 1
            lut_vars[n] = nlutvars

        nr_sim_and_lut_vars = var_idx - 1
        if verbosity > 0:
            print('Created {} SIM and LUT vars'.format(nr_sim_and_lut_vars))

        # Create selection variables. The variable svar_map is a
        # 2-dimensional map from nodesxfanin-size to variablesxnodes.
        # I.e. svar_map[n][1] maps n to a list of 2-tuples which
        # contain all possible options for when node n has 1 fanin. If
        # l = svar_map[n][1] and t = l[0], then t[0] is the variable
        # encoding the situation in which n has fanin t[1].
        svar_map = {} 
        pi_fanin_options = [self.nodes[x] for x in range(self.nr_pis)]

        for n in self.nodes:
            if n.is_pi:
                continue
            svar_map[n] = {}
            fanin_options = None
            if (self.border_io and n.is_border_node) or not self.border_io:
                fanin_options = pi_fanin_options + n.virtual_fanin
            else:
                fanin_options = n.virtual_fanin
            for k in range(min_fanin, max_fanin+1):
                svar_map[n][k] = []
                # Select all possible k-tuples from the fanin options.
                for t in itertools.combinations(fanin_options, k):
                    vt = [0] * (k + 1)
                    vt[0] = var_idx
                    for kp in range(k):
                        vt[kp+1] = t[kp]
                    svar_map[n][k].append(vt)
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
                    rev_var_map[var_idx] = 'houtvars[({},{}] = {}'.format(x, y, var_idx)
                    var_idx += 1
            out_vars[h] = houtvars

        nr_out_vars = var_idx - 1 - nr_sim_and_lut_vars - nr_selection_vars
        if verbosity > 0:
            print('Created {} output vars'.format(nr_out_vars))

        # Create gate-type variables
        gate_type_vars = {}
        for n in self.nodes:
            if n.is_pi:
                continue
            ngatetypevars = [0] * nr_gate_types
            for i in range(nr_gate_types):
                ngatetypevars[i] = var_idx
                rev_var_map[var_idx] = 'ngatetypevars[{}] = {}'.format(i, var_idx)
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
                rev_var_map[var_idx] = 'path_vars[{}][{}] = {}'.format(n, np, var_idx)
                var_idx += 1

        nr_path_vars = var_idx - 1 - nr_sim_and_lut_vars - nr_selection_vars - nr_out_vars - nr_gate_type_vars
        if verbosity > 0:
            print('Created {} path vars'.format(nr_path_vars))

        clauses = []
        # Create simulation propagation constraints
        for n in self.nodes:
            if n.is_pi:
                continue
            for k in range(min_fanin, max_fanin+1):
                fanin_options = svar_map[n][k]
                for option in fanin_options:
                    svar = option[0]
                    for tt_idx in range(nr_sim_vars):
                        permutations = list(itertools.product('01', repeat=(k+1)))
                        for permutation in permutations:
                            function_output = int(permutation[0])
                            const_vals = []
                            f_idx = 0
                            for i in range(k):
                                const_val = int(permutation[i+1])
                                const_vals.append(const_val)
                                f_idx += (const_val << i)

                            clause = [0] * (k + 3)
                            clause[0] = -svar
                            if function_output == 1:
                                clause[1] = -sim_vars[n][tt_idx]
                            else:
                                clause[1] = sim_vars[n][tt_idx]
                            for i in range(len(const_vals)):
                                if const_vals[i] == 1:
                                    # option[i+1] refers to the i-th fanin node
                                    clause[i+2] = -sim_vars[option[i+1]][tt_idx]
                                else:
                                    clause[i+2] = sim_vars[option[i+1]][tt_idx]
                            if function_output == 1:
                                clause[k+2] = lut_vars[n][f_idx]
                            else:
                                clause[k+2] = -lut_vars[n][f_idx]
                            clauses.append(clause)

        # Make sure each node selects exactly one fanin option. Note
        # that the number of fanins it needs to select depends on its
        # gate type: e.g. AND gates need 2 fanins, MAJ gates need
        # 3. If a node selects a k-fanin gate, we disable selection of
        # all options that have more or fewer fanin.
        for n in self.nodes:
            if n.is_pi:
                continue
            ngatetypevars = gate_type_vars[n]
            for i in range(len(enabled_gates)):
                gate_type_var = ngatetypevars[i]
                gate_type = enabled_gates[i]
                fanin_range = self.gate_fanin_range[gate_type]
                kfanin_options = svar_map[n][fanin_range]
                # option[0] is the selection variable whose truth
                # implies that node n selects the fanins nodes in fanin[1:],
                # which should number fanin_range.
                kfanin_vars = [option[0] for option in kfanin_options]
                cnf = CardEnc.equals(lits=kfanin_vars, encoding=EncType.pairwise)
                for clause in cnf.clauses:
                    clauses.append([-gate_type_var] + clause)
                # Selection of this gate type should prevent all fanin
                # options with more or fewer fanin.
                for j in range(len(enabled_gates)):
                    if j == i:
                        continue
                    diff_var = ngatetypevars[j]
                    diff_type = enabled_gates[j]
                    diff_fanin_range = self.gate_fanin_range[diff_type]
                    if diff_fanin_range == fanin_range: # same fanin
                        continue
                    diff_fanin_options = svar_map[n][diff_fanin_range]
                    diff_fanin_vars = [option[0] for option in diff_fanin_options]
                    clauses.append([-gate_type_var, -diff_var])
                    for var in diff_fanin_vars:
                        clauses.append([-gate_type_var, -var])
                        
        # Create cycle-prevention constraints.
        for n in self.nodes:
            if n.is_pi:
                continue
            for k in range(min_fanin, max_fanin+1):
                fanin_options = svar_map[n][k]
                for option in fanin_options:
                    # Suppose node n selects this option, so option[0]
                    # is true.  Then, there is a connection from all
                    # nodes i in option[1:] to n.  Therefore, we must
                    # set the variable encoding path[i][n] to true.
                    # Note that we must filter out PI nodes from this
                    # list.
                    svar = option[0]
                    innodes = [innode for innode in option[1:] if not innode.is_pi]
                    for innode in innodes:
                        clause = [0] * 2
                        clause[0] = -svar
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
            ngatetypevars = gate_type_vars[n]
            fanout_vars = []
            fanout_nodes = []
            for np in self.nodes:
                if np.is_pi or np == n:
                    continue
                for k in range(min_fanin, max_fanin+1):
                    fanin_options = svar_map[np][k]
                    for option in fanin_options:
                        if n in option[1:]:
                            fanout_vars.append(option[0])
                            fanout_nodes.append(np.coords)
            if verbosity > 1:
                print('fanout_vars {}: {}'.format(n.coords, list(zip(fanout_nodes, fanout_vars))))
            # Create cardinality constraints based on gate type.
            for i in range(len(enabled_gates)):
                bound = 1
                if enabled_gates[i] == 'WIRE':                    
                        bound = 3
                gt_var = ngatetypevars[i]
                cnf = CardEnc.atmost(lits=fanout_vars, encoding=EncType.pairwise, bound=bound)
                if verbosity > 1:
                    print('fanout cardinality clauses: {}'.format([[-gt_var] + clause for clause in cnf.clauses]))
                for clause in cnf.clauses:
                    clauses.append([-gt_var] + clause)
                

        # If designated_io is enabled only WIRE elements can have
        # PI/PO fanin/fanout.
        if self.designated_pi:
            assert(self.enable_wire)
            for n in self.nodes:
                if n.is_pi:
                    continue
                if self.border_io and not n.is_border_node:
                    continue
                fanout_vars = []
                fanout_options = []
                for np in self.nodes:
                    if np.is_pi or np == n:
                        continue
                    for k in range(min_fanin, max_fanin+1):
                        fanin_options = svar_map[np][k]
                        for option in fanin_options:
                            if n in option[1:]:
                                fanout_vars.append(option[0])
                                fanout_options.append(np.coords)

                # If enabled, WIRE is the first gate type variable.
                ngatetypevars = gate_type_vars[n]
                wire_type_var = ngatetypevars[0]
                for k in range(min_fanin, max_fanin+1):
                    fanin_options = svar_map[n][k]
                    for option in fanin_options:
                        if any(innode.is_pi for innode in option[1:]):
                            svar = option[0]
                            clauses.append([-svar, wire_type_var])
                            # A designated PI can only have single fanout.
                            cnf = CardEnc.atmost(lits=fanout_vars, encoding=EncType.pairwise)
                            #print('would add: {}'.format(cnf.clauses))
                            for clause in cnf.clauses:
                                clauses.append([-svar] + clause)
    
        if self.designated_po:
            assert(self.enable_wire)
            for n in self.nodes:
                if n.is_pi:
                    continue
                if self.border_io and not n.is_border_node:
                    continue
                fanout_vars = []
                for np in self.nodes:
                    if np.is_pi or np == n:
                        continue
                    for k in range(min_fanin, max_fanin+1):
                        fanin_options = svar_map[np][k]
                        for option in fanin_options:
                            if n in option[1:]:
                                fanout_vars.append(option[0])
                # If enabled, WIRE is the first gate type variable.
                ngatetypevars = gate_type_vars[n]
                wire_type_var = ngatetypevars[0]
                # If one of the POs points to this gate, it has to be
                # a WIRE. Moreover, it cannot have any other fanout.
                for h in range(nr_outputs):
                    houtvars = out_vars[h]
                    houtvar = houtvars[n.coords]
                    clauses.append([-houtvar, wire_type_var])
                    for fanout_var in fanout_vars:
                        clauses.append([-houtvar, -fanout_var])
                        
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
            
            net = logic_network(self.shape, self.nr_pis, self.nr_pos)
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
                fanin_options = svar_map[n][fanin_range]
                nr_valid_fanin_options = 0
                for option in fanin_options:
                    svar = option[0]
                    model_val = model[svar - 1]
                    if model_val > 0:
                        nr_valid_fanin_options += 1
                        fanin_idx = 0
                        for innode in option[1:]:
                            if innode.is_pi:
                                netnode.set_fanin(fanin_idx, net.nodes[innode.coords])
                            else:
                                netnode.set_fanin(fanin_idx, net.node_map[innode.coords])
                            fanin_idx += 1
                assert(nr_valid_fanin_options == 1)

            yield net

