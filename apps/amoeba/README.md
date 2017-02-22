amoeba
======

Amoeba-based superpixel partitioning of multispectral images
---

**About**

Amoeba-based superpixel partitioning of multispectral images into elementary, uniform, connected segments.

**Summary**

We implement an algorithm for scene segmentation that:
* forms a one-to-many partitioning (_"over-segmentation"_) of features in the scene into smaller segments of distinct spectra,
* works essentially like a _kmeans_ based local clustering and aggregating of pixels, and 
* local redundancy in the data, by reducing noise and variability while enforcing connectivity,

so that it can efficiently generate compact, connected, and nearly uniform segments (_"superpixels"_). 

In practice, the approach implemented is based on the estimation of amoeba-like neighborhoods around selected cluster centers that exploit the connections between successive image pixels along geodesic paths in the image. Instead of the pixel grid, it provides a representation of the image into perceptually meaningful entities that can be used as elementary units of any detection, categorization or localization scheme. Indeed, the superpixels capture the spatial/spectral redundancy in images and greatly reduce the complexity of subsequent image processing tasks. They provide convenient primitives from which to compute local image features when objects present in the scene have diverse scales or when they are not known in advance. 

**Description**

Scene segmentation captures the .



**Usage** 

Check the documentation available at: https://gjacopo.github.io/imtools/segmentation/amoebasuperpix.html.

Install Matlab package [`imtools`](https://gjacopo.github.io/imtools/) and run function  [`amoebasuperpix.m`](https://gjacopo.github.io/imtools/segmentation/amoebasuperpix.m).

**<a name="References"></a>References** 

* <a name=“ASSLFS10”></a>R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua and S. Susstrunk (2010): ”SLIC superpixels”, EPFL Technical Report no. 149300.
* Achanta, R., Shaji, A., Smith, K., Lucchi, A., Fua, P., Susstrunk, S., 2012. SLIC
superpixels compared to state-of-the-art superpixel methods. IEEE Transactions
on Pattern Analysis and Machine Intelligence 34 (11), 2274–2282.
* <a name=“FVS09”></a>B. Fulkerson, A. Vedaldi and S. Soatto (2009): ”Class segmentation and object localization with superpixel neighborhoods”, in _Proc. IEEE ICCV_, pp. 670–677, doi:[10.1109/ICCV.2009.5459175](http://dx.doi.org/10.1109/ICCV.2009.5459175).
* <a name=“GS09”></a>Grazzini J. and Soille P. (2009): [**Edge-preserving smoothing using a similarity measure in adaptive geodesic neighbourhoods**](http://www.sciencedirect.com/science/article/pii/S003132030800469X), _Pattern Recognition_, 42(10):2306-2316, doi:[10.1016/j.patcog.2008.11.004](http://dx.doi.org/10.1016/j.patcog.2008.11.004).
* <a name=“GDS10”></a>Grazzini J., Dillard S., and Soille P. (2010): [**Multichannel image regularisation using anisotropic geodesic filtering**](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5596008), in _Proc. IEEE ICPR_, pp. 2664-2667, doi:[10.1109/ICPR.2010.653](http://dx.doi.org/10.1109/ICPR.2010.653).
* <a name=“GS08”></a>Grazzini J.  and Soille P. (2008): [**Adaptive morphological filters using similarities based on geodesic time**](http://www.springerlink.com/content/f6v62233xqkklq72), in _Proc. DGCI_, Lecture Notes in Computer Science, vol. 4992, pp.519-528, doi:[10.1007/978-3-540-79126-3_46](http://dx.doi.org/10.1007/978-3-540-79126-3_46).
* <a name=“Hanbury08”></a> Hanbury A. (2008): How do superpixels affect image segmentation? In: Ruiz-
Shulcloper, J., Kropatsch, W. (Eds.), Progress in Pattern Recognition, Image
Analysis and Applications - Proc. Iberoamerican Congress on Pattern Recognition.
No. 5197 in LNCS. Springer-Verlag, pp. 178–186.
* <a name=“LDM07”></a>R. Lerallut, E. Decenciere, and F. Meyer (2007): ”Image filtering using morphological amoebas”, _Image and Vision Computing_, 25(4):395–404, do:[10.1016/j.imavis.2006.04.018](http://dx.doi.org/10.1016/j.imavis.2006.04.018).
* <a name=“MPWMJ08”></a>A.P. Moore, S. Prince, J. Warrell, U. Mohammed, and G. Jones (2008): ”Superpixel lattices”, in _Proc. IEEE CVPR_, pp. 1–8, doi:[10.1109/CVPR.2008.4587471](http://dx.doi.org/10.1109/CVPR.2008.4587471).
* <a name=“FH04”></a>P. Felzenszwalb and D. Huttenlocher (2004): ”Efficient graph-based image segmentation”, _International Journal of Computer Vision_, 59(2):167-181, doi:[10.1023/B:VISI.0000022288.19776.77](http://dx.doi.org/10.1023/B:VISI.0000022288.19776.77).
