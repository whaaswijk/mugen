Introduction
============

Mugen is a Python library for gate-level SAT-based physical synthesis
of `Field-coupled Nanocomputing
(FCN) <https://www.springer.com/de/book/9783662437216>`_ devices. Its
design is inspired by Marcel Walter's
`fiction <https://github.com/marcelwa/fiction>`_ library. In the
simplest case, Mugen users specify a clocking schemes (encoding
potential gate topologies) and a functional specification (one or more
truth tables encoding circuit functionality). Based on this
information, Mugen returns a logic network satisfying the
specification.

More advanced features include:

- Specifying the gate library. (E.g. including or excluding
  AND/OR/NOT/MAJ gates.)
- Toggling I/O pins.
- Toggling border I/O.

Mugen has been tested under Linux and should generally work on Unix systems.
Windows is not officially supported, although much of the functionality may
still work.
