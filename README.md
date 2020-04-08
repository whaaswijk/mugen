# Mugen

Mugen is a Python library for gate-level SAT-based physical synthesis
for [Field-coupled Nanocomputing
(FCN)](https://www.springer.com/de/book/9783662437216) devices. Its
design is inspired by Marcel Walter's
[fiction](https://github.com/marcelwa/fiction) library. In the
simplest case, Mugen users specify a clocking schemes (encoding
potential gate topologies) and a functional specification (one or more
truth tables encoding circuit functionality). Based on this
information, Mugen returns a logic network satisfying the
specification.

More advanced features include:
- Specifying the gate library. (E.g. including or excluding
  AND/OR/NOT/MAJ gates).
- Support for wire crossings.
- Toggling I/O pins.
- Toggling border I/O.

# Dependencies
## PySAT
Mugen depends on [PySAT](https://github.com/pysathq/pysat). The
easiest way to install PySAT is to simply use pip:
```sh
pip install python-sat
```
## Graphviz (Optional)
This dependency is only required if you want to visualize clocking scheme graphs or logic networks. Installing
the Python package for Graphviz can be done by running:
```sh
pip install graphviz
```
