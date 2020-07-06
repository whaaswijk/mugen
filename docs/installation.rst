.. _installation-label:

Installation
============

The use of Mugen requires the presence of a small number of dependencies.

Numpy/Anaconda
--------------
Mugen depends on `Numpy <https://numpy.org/>`_, and therefore the
easiest way of using Mugen is by installing a scientific Python
distribution first (see
`here <https://www.scipy.org/install.html>`_). The recommended
distribution is `Anaconda <https://www.anaconda.com>`_. Anaconda
will automatically pull in the Numpy dependency as well. 

PySAT
-----
Mugen uses `PySAT <https://github.com/pysathq/pysat>`_ as a SAT
backend. PySAT depends on a **C/C++** compiler as well as the **zlib**
libraries, so just install whatever versions of those work best for
your development environment. The easiest way to install PySAT istself
is to simply use pip:

.. code-block:: sh

    pip install python-sat

wrapt-timeout-decorator 
-----------------------
To support timeouts, Mugen uses
`wrapt-timeout-decorator <https://pypi.org/project/wrapt-timeout-decorator/>`_.
It can be installed using pip:

.. code-block:: sh

    pip install wrapt_timeout_decorator

Glucose::MultiSolvers
---------------------
Mugen can use Glucose::MultiSolvers to perform parallel synthesis. You can
download Glucose 4.x `here <https://www.labri.fr/perso/lsimon/glucose/>`_. In
order for parallel synthesis to work, you must first build the parallel glucose
solver, which can be found under ${GLUCOSE_ROOT}/parallel. You then have to add
${GLUCOSE_ROOT}/parallel to your PATH so that Mugen can find the Glucose binary. 


Graphviz (Optional)
-------------------
This dependency is only required if you want to visualize clocking scheme graphs or logic networks. Installing
the Python package for Graphviz can be done by running:

.. code-block:: sh

    pip install graphviz

The Python interface to Graphviz depends on the **dot** executable
being on your system's path. Install this using the preferred method
for your system. For instance, on a Debian system one might use:

.. code-block:: sh

    sudo apt install graphviz
