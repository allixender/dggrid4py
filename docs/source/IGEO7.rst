IGEO7: A Hierarchically Indexed Hexagonal Equal-Area Discrete Global Grid System
================================================================================

Overview
--------

Hexagonal Discrete Global Grid Systems (DGGS) offer significant advantages for spatial analysis due to their uniform cell shapes and efficient indexing. Among the three central place apertures (3, 4, and 7), aperture 7 subdivisions exhibit very desirable properties, including the preservation of hexagonal symmetry and the formation of unambiguous indexing hierarchies.

Interest in hierarchically indexed aperture 7 hexagonal DGGS has recently increased due to the popularity of the H3 DGGS. But there are currently no open-source equal-area aperture 7 hexagonal DGGS available, that provide similar indexing capabilities like H3.

We present **IGEO7**, a novel pure aperture 7 hexagonal DGGS, and **Z7**, its associated hierarchical integer indexing system. In contrast to H3, where cell sizes vary by up to ±50% across the globe, IGEO7 uses cells of equal area, making it a true equal-area DGGS.

IGEO7 and Z7 are implemented in the open-source software DGGRID. We also present a use case for on-demand suitability modeling to demonstrate a practical application of this new DGGS.


Practical information and OGC compliance
----------------------------------------

The original IGEO7 implementation is available in the DGGRID software since version 8.41 and is a ISEA7H type DGGS with the new Z7 indexing system.

Recent developments in the OGC API DGGS standard draft make reference to refinement ratio 7 DGGRS. In order to have DGGH and ZIRS compliant with the OGC DGGS standard,
we aim to enable a few minor adjustments to IGEO7 through the use of dggrid4py:

-  a rotation of the base icosahedron of 0.5 degrees to align the vertices with water bodies.
- the use of geographiclib to apply authalic conversion to the WGS84 ellipsoid instead of the spherical approximation.

In practice this would mean that IGEO7 as described in the original publication is not the same as the OGC compliant IGEO7 version, though both can be generated through dggrid4py.

We hope to make the required changes available in DGGRID in the future, so that now further confusions can arise. Until then, there might 
be two slightly different IGEO7 implementations in use and implementers shall be explicit. 

API Reference
-------------

The IGEO7 implementation is available through the dggrid4py Python package.

.. seealso::

   :mod:`dggrid4py.igeo7` - IGEO7 module API documentation

Publication
-----------

For more details, see the published paper:

.. seealso::
   
   Kmoch, A., Sahr, K., Chan, W. T., and Uuemaa, E. (2025). IGEO7: A new hierarchically indexed hexagonal equal-area discrete global grid system. *AGILE GIScience Series*, 6, 32. https://doi.org/10.5194/agile-giss-6-32-2025
   
   Full text available at: https://agile-giss.copernicus.org/articles/6/32/2025/


If you use IGEO7 in your research, please cite:

.. code-block:: bibtex

   @article{kmoch2025igeo7,
     author = {Kmoch, A. and Sahr, K. and Chan, W. T. and Uuemaa, E.},
     title = {IGEO7: A new hierarchically indexed hexagonal equal-area discrete global grid system},
     journal = {AGILE GIScience Series},
     volume = {6},
     pages = {32},
     year = {2025},
     doi = {10.5194/agile-giss-6-32-2025},
     url = {https://doi.org/10.5194/agile-giss-6-32-2025}
   }
