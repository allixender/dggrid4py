dggrid4py - a Python library to run highlevel functions of DGGRID
=================================================================

|DOI|

|Population Gridded|

GNU AFFERO GENERAL PUBLIC LICENSE

`DGGRID <https://www.discreteglobalgrids.org/software/>`__ is a free
software program for creating and manipulating Discrete Global Grids
created and maintained by Kevin Sahr. DGGRID version 7.0 was released in
September, 2019.

-  `DGGRID Version 7.5 on GitHub <https://github.com/sahrk/DGGRID>`__
-  `DGGRID User
   Manual <https://raw.githubusercontent.com/sahrk/DGGRID/b61b1d93553c76fc0060f22ed4f903fa2e6789ec/dggridManualV75.pdf>`__

Contents
--------

.. toctree::

   usage
   api


Related work:
-------------

Originally insprired by
`dggridR <https://github.com/r-barnes/dggridR>`__, Richard Barnes’ R
interface to DGGRID. However, dggridR is directly linked via Rcpp to
DGGRID and calls native C/C++ functions.

After some unsuccessful trials with ctypes, cython, CFFI, pybind11 or
cppyy (rather due to lack of experience) I found
`am2222/pydggrid <https://github.com/am2222/pydggrid>`__ (`on
PyPI <https://pypi.org/project/pydggrid/>`__) which made apparently some
initial scaffolding for the transform operation with
`pybind11 <https://pybind11.readthedocs.io/en/master/>`__ including some
sophisticated conda packaging for Windows. This might be worth following
up. Interestingly, its todos include “Adding GDAL export Geometry
Support” and “Support GridGeneration using DGGRID” which this dggrid4py
module supports with integration of GeoPandas.

Bundling for different operating systems
----------------------------------------

Having to compile DGGRID for Windows can be a bit challenging. We are
working on a conda package.

greater context DGGS in Earth Sciences and GIS
----------------------------------------------

Some reading to be excited about:
`discourse.pangeo.io <https://discourse.pangeo.io/t/discrete-global-grid-systems-dggs-use-with-pangeo/2274>`__


Check out the :doc:`usage` section for further information, including
how to :ref:`installation` the project.

.. note::

   This project is under active development.


.. |DOI| image:: docs/source/zenodo.svg
   :target: https://zenodo.org/badge/latestdoi/295495597
.. |Population Gridded| image:: docs/source/day-04-hexa.png
   :target: https://twitter.com/allixender/status/1324055326111485959

