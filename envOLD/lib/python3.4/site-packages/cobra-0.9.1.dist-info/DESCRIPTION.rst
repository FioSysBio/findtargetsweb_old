cobrapy - constraint-based reconstruction and analysis in python
================================================================

|Build Status| |Coverage Status| |Build status| |PyPI| |Gitter| |Waffle|

What is cobrapy?
~~~~~~~~~~~~~~~~

COBRA methods are widely used for genome-scale modeling of metabolic
networks in both prokaryotes and eukaryotes. **cobrapy** is a
constraint-based modeling package that is designed to accommodate the
biological complexity of the next generation of COBRA models and
provides access to commonly used COBRA methods, such as flux balance
analysis, flux variability analysis, and gene deletion analyses.

Our aim with cobrapy is to provide useful, efficient infrastructure
for:

- creating and managing metabolic models
- accessing popular solvers
- analyzing models with methods such as FVA, FBA, pFBA, MOMA etc.
- inspecting models and drawing conclusions on gene essentiality,
  testing consequences of knock-outs etc.

Our goal with cobrapy is for it to be useful on its own, and for it to
be the natural choice of infrastructure for developers that want to
build new COBRA related python packages for e.g. visualization,
strain-design and data driven analysis. By re-using the same classes
and design principles, we can make new methods both easier to
implement and easier to use, thereby bringing the power of COBRA to
more researchers.

The documentation is browseable online at
`readthedocs <https://cobrapy.readthedocs.org/en/stable/>`_ and can
also be
`downloaded <https://readthedocs.org/projects/cobrapy/downloads/>`_.

Please use the `Google
Group <http://groups.google.com/group/cobra-pie>`_ for help. By writing a well formulated question, with sufficient
detail, you are much more likely to quickly receive a good answer! Please refer to these `StackOverflow
guidelines <https://stackoverflow.com/help/how-to-ask>`_ on how to ask questions.
Alternatively, you can use `gitter.im <https://gitter.im/opencobra/cobrapy>`_ for quick questions
and discussions about cobrapy (faster response times). Please keep in mind that answers are provided on a volunteer
basis.

More information about opencobra is available at the
`website <http://opencobra.github.io/>`_.

If you use cobrapy in a scientific publication, please cite
`doi:10.1186/1752-0509-7-74 <http://dx.doi.org/doi:10.1186/1752-0509-7-74>`_

Installation
~~~~~~~~~~~~

Use pip to install cobrapy from
`PyPI <https://pypi.python.org/pypi/cameo>`_ (we recommend doing this
inside a `virtual
environment <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_)::

    pip install cobra

In case you downloaded the source code, run::

    pip install -e .

In the ``cobrapy`` directory. For further information, please follow
the `detailed installation instructions <INSTALL.rst>`_.

Contributing
~~~~~~~~~~~~

Contributions are always welcome! Please read the `contributions
guideline <.github/CONTRIBUTING.rst>`_ to get started.

License
-------

The cobrapy source is released under both the GPL and LGPL licenses. You
may choose which license you choose to use the software under. However,
please note that binary packages which include GLPK (such as the binary
wheels distributed at https://pypi.python.org/pypi/cobra) will be bound
by its license as well.

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License or the Lesser GNU
General Public License as published by the Free Software Foundation,
either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

.. |Build Status| image:: https://travis-ci.org/opencobra/cobrapy.svg?branch=master
   :target: https://travis-ci.org/opencobra/cobrapy
.. |Coverage Status| image:: https://codecov.io/github/opencobra/cobrapy/coverage.svg?branch=master
   :target: https://codecov.io/github/opencobra/cobrapy
.. |Build status| image:: https://ci.appveyor.com/api/projects/status/2o549lhjyukke8nd/branch/master?svg=true
   :target: https://ci.appveyor.com/project/hredestig/cobrapy/branch/master
.. |PyPI| image:: https://img.shields.io/pypi/v/cobra.svg
   :target: https://pypi.python.org/pypi/cobra
.. |Gitter| image:: https://badges.gitter.im/opencobra/cobrapy.svg
   :target: https://gitter.im/opencobra/cobrapy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
.. |Waffle| image:: https://badge.waffle.io/opencobra/cobrapy.png?label=ready&title=Ready
   :target: https://waffle.io/opencobra/cobrapy
   :alt: 'Stories in Ready'


Installation of cobrapy
=======================

For installation help, please use the `Google
Group <http://groups.google.com/group/cobra-pie>`_. For usage
instructions, please see the
`documentation <https://cobrapy.readthedocs.org/en/latest/>`_.

All releases require Python 2.7+ or 3.4+ to be installed before
proceeding. Mac OS X (10.7+) and Ubuntu ship with Python. Windows users
without python can download and install python from the `python
website <https://www.python.org/ftp/python/2.7.9/python-2.7.9.amd64.msi>`_.
Please note that though Anaconda and other python distributions may work
with cobrapy, they are not explicitly supported (yet!).

Stable version installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

cobrapy can be installed with any recent installation of pip.
Instructions for several operating systems are below:

Mac OS X or Linux
-----------------

1. `install
   pip <http://pip.readthedocs.org/en/latest/installing.html>`_.
2. In a terminal, run ``sudo pip install cobra``

We highly recommend updating ``pip`` beforehand (``pip install pip --upgrade``).

Microsoft Windows
-----------------

The preferred installation method on Windows is also to use pip. The
latest Windows installers for Python 2.7 and 3.4 include pip, so if you
use those you will already have pip.

1. In a terminal, run ``C:\Python27\Scripts\pip.exe install cobra`` (you
   may need to adjust the path accordingly).

To install without pip, you will need to download and use the
appropriate installer for your version of python from the `python
package index <https://pypi.python.org/pypi/cobra/>`_.

Installation for development
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get the `detailed contribution instructions <CONTRIBUTING.rst>`_ for contributing to cobrapy.

Installation of optional dependencies
=====================================

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

On windows, these can downloaded from [this site]
(http://www.lfd.uci.edu/~gohlke/pythonlibs/). On Mac/Linux, they can be
installed using pip, or from the OS package manager (e.g brew, apt,
yum).

1. `libsbml <http://sbml.org>`_ >= 5.10 to read/write SBML level 2
   files

   -  `Windows libsbml installer <http://www.lfd.uci.edu/~gohlke/pythonlibs/#libsbml>`_
   -  Use ``sudo pip install python-libsbml`` on Mac/Linux

2. `lxml <http://lxml.de/>`_ to speed up read/write of SBML level 3 files.

3. `scipy <http://scipy.org>`_ >= 0.11 for MOMA and saving to \*.mat files.
   -  `Windows scipy installer <http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy>`_

4. `pytest <http://docs.pytest.org/en/latest/>`_ and `pytest-benchmark
   <http://pytest-benchmark.readthedocs.io/en/latest/>`_ are required
   for testing

You can install all packages directly by

.. code:: shell

    pip install "cobra[all]"

Solvers
~~~~~~~~~~~~~

cobrapy uses `optlang <http://optlang.readthedocs.io>`_ to interface
the mathematical solvers used to optimize the created COBRA models,
which at the time of writing

-  ILOG/CPLEX (available with
   `Academic <https://www.ibm.com/developerworks/university/academicinitiative/>`_
   and
   `Commercial <http://www.ibm.com/software/integration/optimization/cplex-optimizer/>`_
   licenses).
-  `gurobi <http://gurobi.com>`_
-  `glpk <http://www.gnu.org/software/glpk/>`_

Testing your installation
=========================

While it is not a hard requirement for using cobrapy, you need pytest and
pytest-benchmark to run its tests. First do

.. code:: shell

    pip install pytest pytest-benchmark

or to install cobrapy directly with the test dependencies

.. code:: shell

   pip install "cobra[test]"

Then start python and type the following into the Python shell

.. code:: python

    from cobra.test import test_all
    test_all()

You should see some skipped tests and expected failures, and the
function should return ``0``. If you see a value other than ``0`` please file
an `issue report <CONTRIBUTING.rst>`_.


