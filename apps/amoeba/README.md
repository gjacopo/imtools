amoeba
======

Amoeba-based superpixel partitioning of multispectral images
---

**Summary**

We implement an algorithm for scene segmentation that:
* forms a one-to-many partitioning (_"over-segmentation"_) of features in the scene into smaller segments of distinct spectra,
* works essentially like a _kmeans_ based local clustering and aggregating of pixels, and 
* local redundancy in the data, by reducing noise and variability while enforcing connectivity,

so that it can efficiently generate compact, connected, nearly uniform and erceptually meaningful segments (_"superpixels"_). 

In practice, the implementation is based on the estimation of amoeba-like neighborhoods around selected cluster centers that exploit the connections between successive image pixels along geodesic paths in the image. The resulting superpixels capture the spatial/spectral redundancy and greatly reduce the complexity of subsequent image processing tasks. They provide convenient primitives from which to compute local image features when objects present in the scene have diverse scales or when they are not known in advance. 

**Description**


**Usage** 

Check the documentation available at: https://gjacopo.github.io/imtools/segmentation/amoebasuperpix.html.

Install Matlab package [`imtools`](https://gjacopo.github.io/imtools/) and run function  [`amoebasuperpix.m`](https://gjacopo.github.io/imtools/segmentation/amoebasuperpix.m).

**<a name="References"></a>References** 

* <a name=“ASSLFS10”></a>R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua and S. Susstrunk (2010): ”SLIC superpixels”, EPFL Technical Report no. 149300.
* <a name=“GS09”></a>Grazzini J. and Soille P. (2009): [**Edge-preserving smoothing using a similarity measure in adaptive geodesic neighbourhoods**](http://www.sciencedirect.com/science/article/pii/S003132030800469X), _Pattern Recognition_, 42(10):2306-2316, doi:[10.1016/j.patcog.2008.11.004](http://dx.doi.org/10.1016/j.patcog.2008.11.004).
* <a name=“GDS10”></a>Grazzini J., Dillard S., and Soille P. (2010): [**Multichannel image regularisation using anisotropic geodesic filtering**](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5596008), in _Proc. IEEE ICPR_, pp. 2664-2667, doi:[10.1109/ICPR.2010.653](http://dx.doi.org/10.1109/ICPR.2010.653).
* <a name=“GS08”></a>Grazzini J.  and Soille P. (2008): [**Adaptive morphological filters using similarities based on geodesic time**](http://www.springerlink.com/content/f6v62233xqkklq72), in _Proc. DGCI_, Lecture Notes in Computer Science, vol. 4992, pp.519-528, doi:[10.1007/978-3-540-79126-3_46](http://dx.doi.org/10.1007/978-3-540-79126-3_46).
* <a name=“LDM07”></a>R. Lerallut, E. Decenciere, and F. Meyer (2007): ”Image filtering using morphological amoebas”, _Image and Vision Computing_, 25(4):395–404, do:[10.1016/j.imavis.2006.04.018](http://dx.doi.org/10.1016/j.imavis.2006.04.018).
* <a name=“MPWMJ08”></a>A.P. Moore, S. Prince, J. Warrell, U. Mohammed, and G. Jones (2008): ”Superpixel lattices”, in _Proc. IEEE CVPR_, pp. 1–8, doi:[10.1109/CVPR.2008.4587471](http://dx.doi.org/10.1109/CVPR.2008.4587471).
