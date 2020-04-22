# Mugen

Mugen is a Python library for gate-level SAT-based physical synthesis
of [Field-coupled Nanocomputing
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
- Toggling I/O pins.
- Toggling border I/O.

# Dependencies

## Numpy/Anaconda
Mugen depends on [Numpy](https://numpy.org/), and therefore the
easiest way of using Mugen is by installing a scientific Python
distribution first (see
[here](https://www.scipy.org/install.html)). The recommended
distribution is [Anaconda](https://www.anaconda.com). Anaconda
will automatically pull in the Numpy as well. 

## PySAT

Mugen uses [PySAT](https://github.com/pysathq/pysat) as a SAT
backend. PySAT depends on a **C/C++** compiler as well as the **zlib**
libraries, so just install whatever versions of those work best for
your development environment. The easiest way to install PySAT istself
is to simply use pip:

```sh
pip install python-sat
```

## Graphviz (Optional)

This dependency is only required if you want to visualize clocking scheme graphs or logic networks. Installing
the Python package for Graphviz can be done by running:
```sh
pip install graphviz
```

The Python interface to Graphviz depends on the **dot** executable
being on your system's path. Install this using the preferred method
for your system. For instance, on a Debian system one might use:

```sh
sudo apt install graphviz
```
