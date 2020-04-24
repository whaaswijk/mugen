import unittest
from graph import *

class graph_tests(unittest.TestCase):

    def test_graph_construction(self):
        '''
        Build a simple graph and render it to PNG to make sure the
        construction is correct. This example uses a 3x3 slice from
        the USE topology.
        '''
        # First, create the graph and give it's shape
        g = scheme_graph(shape=(3,3))

        # Next, create virtual connections between nodes/tiles.
        # These are used to specify potential connections that the
        # SAT solver should look at.
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (2, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        g.add_virtual_edge((1, 1), (1, 2))
        g.add_virtual_edge((2, 1), (2, 0))
        g.add_virtual_edge((2, 1), (1, 1))
        g.add_virtual_edge((0, 2), (0, 1))
        g.add_virtual_edge((0, 2), (1, 2))
        g.add_virtual_edge((1, 2), (2, 2))
        g.add_virtual_edge((2, 2), (2, 1))
        
        g.to_png('test')

    def test_3x3_cycles(self):
        '''
        Builds a graph with a 3x3 USE topology and computes its cycles.
        '''
        g = scheme_graph(shape=(3,3))
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (2, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        g.add_virtual_edge((1, 1), (1, 2))
        g.add_virtual_edge((2, 1), (2, 0))
        g.add_virtual_edge((2, 1), (1, 1))
        g.add_virtual_edge((0, 2), (0, 1))
        g.add_virtual_edge((0, 2), (1, 2))
        g.add_virtual_edge((1, 2), (2, 2))
        g.add_virtual_edge((2, 2), (2, 1))

        cycles = g.find_cycles()
        nr_cycles = len(cycles)
#        for cycle in cycles:
#            print(cycle)
        self.assertEqual(nr_cycles, 8)

    def test_4x4_cycles(self):
        '''
        Builds a graph with a 4x4 USE topology and computes its cycles.
        '''
        g = scheme_graph(shape=(4,4))
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (2, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((2, 0), (3, 0))
        g.add_virtual_edge((3, 0), (3, 1))

        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        g.add_virtual_edge((1, 1), (1, 2))
        g.add_virtual_edge((2, 1), (2, 0))
        g.add_virtual_edge((2, 1), (1, 1))
        g.add_virtual_edge((3, 1), (2, 1))
        g.add_virtual_edge((3, 1), (3, 2))

        g.add_virtual_edge((0, 2), (0, 1))
        g.add_virtual_edge((0, 2), (1, 2))
        g.add_virtual_edge((1, 2), (2, 2))
        g.add_virtual_edge((1, 2), (1, 3))
        g.add_virtual_edge((2, 2), (2, 1))
        g.add_virtual_edge((2, 2), (3, 2))
        g.add_virtual_edge((3, 2), (3, 3))

        g.add_virtual_edge((0, 3), (0, 2))
        g.add_virtual_edge((1, 3), (0, 3))
        g.add_virtual_edge((2, 3), (1, 3))
        g.add_virtual_edge((2, 3), (2, 2))
        g.add_virtual_edge((3, 3), (2, 3))

        cycles = g.find_cycles()
        nr_cycles = len(cycles)
#        for cycle in cycles:
#            print(cycle)
        self.assertEqual(nr_cycles, 176)

        
        
class synth_tests(unittest.TestCase):

    def test_and_synthesis(self):
        '''
        Creates some simple clock scheme graph toplogies and tries
        to synthesize an AND function.
        '''
        # Create a very simple graph for a 1x1 grid with 2 PIs.
        # Disable all gate types except for and. Try to synthesize an
        # AND function and count the number of solutions. There should
        # be a unique solution.
        g = scheme_graph(shape=(1,1))
        g.enable_or = False
        g.enable_not = False
        g.enable_maj = False
        functions = [[0,0,0,1]]

        models_found = 0
        for net in g.synthesize(functions):
            self.assertTrue(g.satisfies_spec(net))
            models_found += 1
        self.assertEqual(models_found, 1)

        # Do the same thing but now say there are 3 PIs. There should
        # still only be one model.
        g = scheme_graph(shape=(1,1))
        g.enable_or = False
        g.enable_not = False
        g.enable_maj = False
        functions = [[0,0,0,1,0,0,0,1]]
        models_found = 0
        for net in g.synthesize(functions):
            self.assertTrue(g.satisfies_spec(net))
            models_found += 1
        self.assertEqual(models_found, 1)

        # Again do the same thing but now allow all gate types. The
        # number of solutions shouldn't change, because MAJ and OR
        # gates do not help us here.
        g = scheme_graph(shape=(1,1))
        functions = [[0,0,0,1,0,0,0,1]]
        models_found = 0
        for net in g.synthesize(functions):
            self.assertTrue(g.satisfies_spec(net))
            models_found += 1
        self.assertEqual(models_found, 1)

    def test_maj_synthesis(self):
        '''
        Test synthesis of a single majority gate with 3 inputs.  There
        should be a unique solution.
        '''
        g = scheme_graph(shape=(1,1))
        functions = [[0,0,0,1,0,1,1,1]]
        models_found = 0
        for net in g.synthesize(functions):
            self.assertTrue(g.satisfies_spec(net))
            models_found += 1
        self.assertEqual(models_found, 1)

    def test_maj_fail(self):
        # If we disable MAJ gates, we shouldn't be able to synthesize
        # a majority function with just one gate.
        g = scheme_graph(shape=(1,1))
        g.enable_maj = False
        functions = [[0,0,0,1,0,1,1,1]]
        models_found = 0
        for net in g.synthesize(functions):
            models_found += 1
        self.assertEqual(models_found, 0)


    def test_and_or_maj(self):
        '''
        Build a circuit that computes a majority function without
        using majority gates.  This is also a good test to make sure
        we don't allow cycles. This should find solutions for 3x3
        USE topologies, but not for 2x2 USE topologies.
        '''
        g = scheme_graph(shape=(3,3))
        g.enable_maj = False
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (2, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        g.add_virtual_edge((1, 1), (1, 2))
        g.add_virtual_edge((2, 1), (2, 0))
        g.add_virtual_edge((2, 1), (1, 1))
        g.add_virtual_edge((0, 2), (0, 1))
        g.add_virtual_edge((0, 2), (1, 2))
        g.add_virtual_edge((1, 2), (2, 2))
        g.add_virtual_edge((2, 2), (2, 1))
        functions = [[0,0,0,1,0,1,1,1]]
        models_found = 0
        for net in g.synthesize(functions):
            self.assertTrue(g.satisfies_spec(net))
            models_found += 1
            if models_found > 1000:
                break
        self.assertTrue(models_found > 0)

        g = scheme_graph(shape=(2,2))
        g.enable_maj = False
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        functions = [[0,0,0,1,0,1,1,1]]
        models_found = 0
        cycles = g.find_cycles()
        for net in g.synthesize(functions):
            models_found += 1
            break;
        self.assertEqual(models_found, 0)

    def test_and_or_mux(self):
        '''
        Build a circuit that computes a MUX function without
        using majority gates. Using a 3x3 USE topology.
        '''
        g = scheme_graph(shape=(3,3))
        g.enable_maj = False
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (2, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        g.add_virtual_edge((1, 1), (1, 2))
        g.add_virtual_edge((2, 1), (2, 0))
        g.add_virtual_edge((2, 1), (1, 1))
        g.add_virtual_edge((0, 2), (0, 1))
        g.add_virtual_edge((0, 2), (1, 2))
        g.add_virtual_edge((1, 2), (2, 2))
        g.add_virtual_edge((2, 2), (2, 1))
        functions = [[0,0,1,1,0,1,0,1]]
        models_found = 0
        for net in g.synthesize(functions): #, verbosity=2):
            self.assertTrue(g.satisfies_spec(net))
            models_found += 1
            if models_found > 1000:
                break
        self.assertTrue(models_found > 0)

    def test_disable_wire(self):
        '''
        Build a circuit that computes a MUX function without using
        majority gates or WIRE elements. Using a 3x3 USE topology.
        '''
        g = scheme_graph(shape=(3,3))
        g.enable_maj = False
        g.enable_wire = False
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (2, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        g.add_virtual_edge((1, 1), (1, 2))
        g.add_virtual_edge((2, 1), (2, 0))
        g.add_virtual_edge((2, 1), (1, 1))
        g.add_virtual_edge((0, 2), (0, 1))
        g.add_virtual_edge((0, 2), (1, 2))
        g.add_virtual_edge((1, 2), (2, 2))
        g.add_virtual_edge((2, 2), (2, 1))
        functions = [[0,0,1,1,0,1,0,1]]
        models_found = 0
        for net in g.synthesize(functions): #, verbosity=2):
            self.assertTrue(g.satisfies_spec(net))
            models_found += 1
            if models_found >= 10:
                break
        self.assertTrue(models_found > 0)

    def test_border_io(self):
        '''
        Build a circuit that computes a MUX function without
        using majority gates. Using a 3x3 USE topology.
        Verifies that only border nodes have I/O pins.
        '''
        g = scheme_graph(shape=(3,3))
        g.enable_maj = False
        g.border_io = True
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (2, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        g.add_virtual_edge((1, 1), (1, 2))
        g.add_virtual_edge((2, 1), (2, 0))
        g.add_virtual_edge((2, 1), (1, 1))
        g.add_virtual_edge((0, 2), (0, 1))
        g.add_virtual_edge((0, 2), (1, 2))
        g.add_virtual_edge((1, 2), (2, 2))
        g.add_virtual_edge((2, 2), (2, 1))
        functions = [[0,0,1,1,0,1,0,1]]
        models_found = 0
        for net in g.synthesize(functions): #, verbosity=2):
            self.assertTrue(g.satisfies_spec(net))
            models_found += 1
            if models_found >= 10000:
                net.to_png('border-io-test')
                break
#        print('Found {} border I/O models'.format(models_found))
        self.assertTrue(models_found > 0)

    def test_designated_pi(self):
        '''
        Builds a circuit that computes a 3-input MUX without using
        majority gates and while enabling "designated_pi". Using a 5x3
        USE topology. Verifies that only WIRE elements have PI fanin.
        '''
        g = scheme_graph(shape=(5,3))
        g.enable_maj = False
        g.designated_pi = True
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (2, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((2, 0), (3, 0))
        g.add_virtual_edge((3, 0), (3, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        g.add_virtual_edge((1, 1), (1, 2))
        g.add_virtual_edge((2, 1), (2, 0))
        g.add_virtual_edge((2, 1), (1, 1))
        g.add_virtual_edge((3, 1), (2, 1))
        g.add_virtual_edge((3, 1), (3, 2))
        g.add_virtual_edge((0, 2), (0, 1))
        g.add_virtual_edge((0, 2), (1, 2))
        g.add_virtual_edge((1, 2), (2, 2))
        g.add_virtual_edge((2, 2), (2, 1))
        g.add_virtual_edge((2, 2), (3, 2))

        g.add_virtual_edge((3, 0), (4, 0))
        g.add_virtual_edge((4, 1), (4, 0))
        g.add_virtual_edge((4, 1), (3, 1))
        g.add_virtual_edge((4, 2), (4, 1))
        g.add_virtual_edge((3, 2), (4, 2))
        
        functions = [[0,0,1,1,0,1,0,1]]
        models_found = 0
        for net in g.synthesize(functions):
            self.assertTrue(g.satisfies_spec(net))
            models_found += 1
            if models_found >= 1000:
                break
#        print('Found {} designated PI models'.format(models_found))            
        self.assertTrue(models_found > 0)

    def test_designated_po(self):
        '''
        Builds a circuit that computes a 3-input MUX without using
        majority gates and while enabling "designated_po". Using a 3x3
        USE topology.  Verifies that only WIRE elements have PO
        fanout.
        '''
        g = scheme_graph(shape=(3,3))
        g.enable_maj = False
        g.designated_po = True
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (2, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        g.add_virtual_edge((1, 1), (1, 2))
        g.add_virtual_edge((2, 1), (2, 0))
        g.add_virtual_edge((2, 1), (1, 1))
        g.add_virtual_edge((0, 2), (0, 1))
        g.add_virtual_edge((0, 2), (1, 2))
        g.add_virtual_edge((1, 2), (2, 2))
        g.add_virtual_edge((2, 2), (2, 1))
        functions = [[0,0,1,1,0,1,0,1]]
        models_found = 0
        for net in g.synthesize(functions): #, verbosity=2):
            self.assertTrue(g.satisfies_spec(net))
            models_found += 1
            if models_found >= 1000:
                break
#        print('Found {} designated PO models'.format(models_found))
        self.assertTrue(models_found > 0)

    def test_restricted_io(self):
        '''
        Builds a circuit that computes a 3-input MUX without
        using majority gates and while enabling "designated_pi",
        "designated_po", and "border_io". Using a 5x3 USE topology.
        Verifies that only WIRE elements on the border have PI/PO
        fanin/fanout.
        '''
        g = scheme_graph(shape=(5,3))
        g.enable_maj = False
        g.designated_pi = True
        g.designated_po = True
        g.border_io = True

        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (2, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((2, 0), (3, 0))
        g.add_virtual_edge((3, 0), (3, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        g.add_virtual_edge((1, 1), (1, 2))
        g.add_virtual_edge((2, 1), (2, 0))
        g.add_virtual_edge((2, 1), (1, 1))
        g.add_virtual_edge((3, 1), (2, 1))
        g.add_virtual_edge((3, 1), (3, 2))
        g.add_virtual_edge((0, 2), (0, 1))
        g.add_virtual_edge((0, 2), (1, 2))
        g.add_virtual_edge((1, 2), (2, 2))
        g.add_virtual_edge((2, 2), (2, 1))
        g.add_virtual_edge((2, 2), (3, 2))
        g.add_virtual_edge((3, 0), (4, 0))
        g.add_virtual_edge((4, 1), (4, 0))
        g.add_virtual_edge((4, 1), (3, 1))
        g.add_virtual_edge((4, 2), (4, 1))
        g.add_virtual_edge((3, 2), (4, 2))

        functions = [[0,0,1,1,0,1,0,1]]
        models_found = 0
        for net in g.synthesize(functions):
            self.assertTrue(g.satisfies_spec(net))
            models_found += 1
            if models_found >= 1000:
                net.to_png('restricted-io')
                break
#        print('Found {} restricted I/O models'.format(models_found))
        self.assertTrue(models_found > 0)

    def test_no_multi_edges(self):
        '''
        Builds a circuit that should not contain multi-edges. Using a 3x3 USE topology.
        '''
        g = scheme_graph(shape=(3,3))
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (2, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        g.add_virtual_edge((1, 1), (1, 2))
        g.add_virtual_edge((2, 1), (2, 0))
        g.add_virtual_edge((2, 1), (1, 1))
        g.add_virtual_edge((0, 2), (0, 1))
        g.add_virtual_edge((0, 2), (1, 2))
        g.add_virtual_edge((1, 2), (2, 2))
        g.add_virtual_edge((2, 2), (2, 1))
        functions = [[0,0,0,1,0,0,0,1,0,0,0,1,0,1,1,1,0,0,0,1,0,1,1,1,0,1,1,1,0,1,1,1]]
        models_found = 0
        for net in g.synthesize(functions): #, verbosity=2):
            models_found += 1
            if models_found >= 1000:
                net.to_png('multi-edges')
                break
        self.assertTrue(models_found > 0) 
        
    def test_XOR2(self):
        '''
        Builds a circuit that computes the XOR-2 function. Using a 3x3 USE topology.
        '''
        g = scheme_graph(shape=(3,3))
        g.enable_maj = False
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (2, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        g.add_virtual_edge((1, 1), (1, 2))
        g.add_virtual_edge((2, 1), (2, 0))
        g.add_virtual_edge((2, 1), (1, 1))
        g.add_virtual_edge((0, 2), (0, 1))
        g.add_virtual_edge((0, 2), (1, 2))
        g.add_virtual_edge((1, 2), (2, 2))
        g.add_virtual_edge((2, 2), (2, 1))
        functions = [[0,1,1,0]]
        models_found = 0
        for net in g.synthesize(functions): #, verbosity=2):
            models_found += 1
            if models_found >= 1000:
                net.to_png('XOR-2')
                break
        self.assertTrue(models_found > 0)
        
    def test_XNOR2(self):
        '''
        Builds a circuit that computes the XOR-2 function. Using a 3x3 USE topology.
        '''
        g = scheme_graph(shape=(3,3))
        g.enable_maj = False
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (2, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        g.add_virtual_edge((1, 1), (1, 2))
        g.add_virtual_edge((2, 1), (2, 0))
        g.add_virtual_edge((2, 1), (1, 1))
        g.add_virtual_edge((0, 2), (0, 1))
        g.add_virtual_edge((0, 2), (1, 2))
        g.add_virtual_edge((1, 2), (2, 2))
        g.add_virtual_edge((2, 2), (2, 1))
        functions = [[1,0,0,1]]
        models_found = 0
        for net in g.synthesize(functions): #, verbosity=2):
            models_found += 1
            if models_found >= 1000:
                net.to_png('XNOR-2')
                break
        self.assertTrue(models_found > 0)
        
    def test_MUX21(self):
        '''
        Builds a circuit that computes the 2:1 MUX function. Using a 3x3 USE topology.
        Designated I/O pins are required at the layout borders.
        '''
        g = scheme_graph(shape=(3,3), border_io=True, designated_pi=True, designated_po=True)
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (2, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        g.add_virtual_edge((1, 1), (1, 2))
        g.add_virtual_edge((2, 1), (2, 0))
        g.add_virtual_edge((2, 1), (1, 1))
        g.add_virtual_edge((0, 2), (0, 1))
        g.add_virtual_edge((0, 2), (1, 2))
        g.add_virtual_edge((1, 2), (2, 2))
        g.add_virtual_edge((2, 2), (2, 1))
        functions = [[0,1,0,1,0,0,1,1]]
        models_found = 0
        for net in g.synthesize(functions): #, verbosity=2):
            models_found += 1
            if models_found >= 1000:
                net.to_png('MUX21')
                break
        self.assertTrue(models_found > 0)

if __name__ == '__main__':
    unittest.main()
