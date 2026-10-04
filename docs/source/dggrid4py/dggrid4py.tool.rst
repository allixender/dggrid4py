dggrid4py.tool
==============

Portable DGGRID binary
----------------------

If DGGRID can neither be installed from conda-forge nor compiled, dggrid4py can download a portable binary from the
`DGGRID_portables <https://github.com/allixender/DGGRID_portables>`_ releases. The binaries are built without GDAL
and are still experimental, see :doc:`../usage` for an example and the points to consider.


.. automodule:: dggrid4py.tool
   :members:

   .. rubric:: Functions

   .. autosummary::
   
      get_portable_executable
      portable_asset_name
      download_executable