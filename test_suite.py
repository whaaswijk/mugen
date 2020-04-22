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

class synth_tests(unittest.TestCase):

    def test_and_synthesis(self):
        '''
        Creates some simple clock scheme graph toplogies and tries
        to synthesize an AND function.
        '''
        # Create a very simple graph for a 1x1 grid with 2 PIs.
        # Disable all gate types except for and. Try to synthesize an
        # AND function and count the number of solutions.
        g = scheme_graph(shape=(1,1))
        g.enable_or = False
        g.enable_not = False
        g.enable_maj = False
        functions = [[0,0,0,1]]

        models_found = 0
        for net in g.synthesize(functions):
            models_found += 1
        self.assertEqual(models_found, 4)

        # Do the same thing but now say there are 3 PIs. The number of
        # solutions will increase because of the extra degree of
        # freedom for the 3rd fanin.
        g = scheme_graph(shape=(1,1))
        g.enable_or = False
        g.enable_not = False
        g.enable_maj = False
        functions = [[0,0,0,1,0,0,0,1]]
        models_found = 0
        for net in g.synthesize(functions):
            models_found += 1
        self.assertEqual(models_found, 6)

        # Again do the same thing but now allow all gate types. The
        # number of solutions shouldn't change, because MAJ and OR
        # gates do not help us here.
        g = scheme_graph(shape=(1,1))
        functions = [[0,0,0,1,0,0,0,1]]
        models_found = 0
        for net in g.synthesize(functions):
            models_found += 1
        self.assertEqual(models_found, 6)

    def test_maj_synthesis(self):
        '''
        Test synthesis of a single majority gate with 3 inputs.  The
        number of solutions should be 3! (=6), since there is a
        single gate and we can pick any PI for any gate fanin.
        '''
        g = scheme_graph(shape=(1,1))
        functions = [[0,0,0,1,0,1,1,1]]
        models_found = 0
        for net in g.synthesize(functions):
            models_found += 1
        self.assertEqual(models_found, 6)

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
            models_found += 1
            break
        self.assertEqual(models_found, 1)

        g = scheme_graph(shape=(2,2))
        g.enable_maj = False
        g.add_virtual_edge((0, 0), (1, 0))
        g.add_virtual_edge((1, 0), (1, 1))
        g.add_virtual_edge((0, 1), (0, 0))
        g.add_virtual_edge((1, 1), (0, 1))
        functions = [[0,0,0,1,0,1,1,1]]
        models_found = 0
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
            models_found += 1
#            print(net)
#            net.to_png('mux-test')
#            input()
            break
        self.assertEqual(models_found, 1)

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
            models_found += 1
            for n in net.nodes:
                self.assertFalse(n.gate_type in ['WIRE', 'MAJ'])
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
            models_found += 1
            if not net.has_border_io():
                net.to_png('border-io-test')
                self.assertTrue(False)
            if models_found > 10000:
                net.to_png('border-io-test')
                break
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
            models_found += 1
            self.assertTrue(net.has_designated_pi())
            if models_found > 1000:
                break
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
            models_found += 1
            self.assertTrue(net.has_designated_po())
            if models_found > 1000:
                break
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
            models_found += 1
            self.assertTrue(net.has_designated_pi() and net.has_designated_po() and net.has_border_io())
            if models_found > 1000:
                break
        self.assertTrue(models_found > 0)

if __name__ == '__main__':
    unittest.main()
