amoeba
======

Amoeba-based superpixel partitioning of multispectral images
---

**About**

Amoeba-based superpixel partitioning of multispectral images into elementary, uniform, connected units.

<table align="center">
    <tr> <td align="left"><i>documentation</i></td> <td align="left">available at: https://gjacopo.github.io/imtools/segmentation/amoebasuperpix.html</td> </tr> 
    <tr> <td align="left"><i>version</i></td> <td align="left">1.0 <i>(non-active development)</i> </td> </tr> 
    <tr> <td align="left"><i>since</i></td> <td align="left">2007</td> </tr> 
    <tr> <td align="left"><i>license</i></td> <td align="left"><a href="https://joinup.ec.europa.eu/sites/default/files/eupl1.1.-licence-en_0.pdfEUPL">EUPL</a>  <i>(cite the source code or any of the references below!)</i> </td> </tr> 
</table>

**Summary**

Segmentation is a difficult task because of the high complexity of images, where complexity refers to the large variety of pictorial representations of objects with the same semantic meaning and also to the extensive amount of available details.
It seems therefore natural, and presumably more efficient, to work with perceptually meaningful entities obtained from low-level grouping processes instead of the pixel representation. In that context, superpixels obtained from conservative over-segmentation are a common pre-processing step for recovering image features. 

We implement a new algorithm that:
* works essentially like a k-means based local clustering of pixels, but 
* enforces connectivity, so that it can efficiently generate compact, connected, and nearly uniform superpixels. 

This approach is based on the estimation of amoeba-like neighborhoods around selected cluster centers that exploit the connections between successive image pixels along geodesic paths in the image. The resulting superpixels capture the spatial/spectral redundancy in images and greatly reduce the complexity of subsequent image processing tasks. They provide convenient primitives from which to compute local image features when objects present in the scene have diverse scales or when they are not known in advance. 

**Description**

The problem of segmenting an image into semantically meaningful units is a fundamental one in image processing. Indeed, visual information extraction is usually regarded as a segmentation issue [[FH04]](FH04). Besides practical issues (running time, user interaction, etc...), the key issue is what objects, if any, correspond to the segmented regions. It is well known that segmentation is an ill-posed problem whose ’correct’ solution is largely dependent on the application, if not completely subjective. Scene segmentation captures the local redundancy in the data, by reducing noise and variability, but aggregating pixels into segments often entails a decision that is unrelated to the final task. The goal is to perform this decision in a conservative way to minimize the risk of merging unrelated pixels from different objects. Over-segmentation is a pragmatic
alternative that forms a one-to-many partitioning of scene features into smaller segments of distinct spectra. Instead of the pixel grid, it provides a representation of the image into perceptually meaningful entities, the so-called superpixels, that can be used as elementary units of any detection, categorization or localization scheme [[FH04]](FH04) [[FVS09]](FVS09) [[MPWMJ08]](MPWMJ08).

**Usage** 

Install Matlab package [`imtools`](https://gjacopo.github.io/imtools/) and run function  [`amoebasuperpix.m`](https://gjacopo.github.io/imtools/segmentation/amoebasuperpix.m).

**<a name="References"></a>References** 

# <a name=“ASSLFS10”></a>R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua and S. Susstrunk (2010): ”SLIC superpixels”, EPFL Technical Report no. 149300.
Achanta, R., Shaji, A., Smith, K., Lucchi, A., Fua, P., Susstrunk, S., 2012. SLIC
superpixels compared to state-of-the-art superpixel methods. IEEE Transactions
on Pattern Analysis and Machine Intelligence 34 (11), 2274–2282.
Bagon, S., Boiman, O., Irani, M., 2008. What is a good image segment? A unified
approach to segment extraction. In: Proc. European Conference on Computer
Vision. Vol. 5305 of Lecture Notes in Computer Science. Springer-Verlag, pp.
30–44.
Bertelli, L., Manjunath, B., 2006. Redundancy in all pairs fast marching method.
In: Proc. IEEE International Conference on Image Processing. pp. 3033–3036.
Coeurjolly, D., Miguet, D., Tougne, L., 2004. 2D and 3D visibility in discrete
geometry: an application to discrete geodesic paths. Pattern Recognition Letters
25 (5), 561–570.
Debayle, J., Pinoli, J., 2006a. General adaptive neighborhood image processing.
Part I: Introduction and theoretical aspects. Journal of Mathematical Imaging
and Vision 25 (2), 245–266.
Debayle, J., Pinoli, J., 2006b. General adaptive neighborhood image processing.
Part II: practical application examples. Journal of Mathematical Imaging and
Vision 25 (2), 267–284.
# <a name=“FVS09”></a>B. Fulkerson, A. Vedaldi and S. Soatto (2009): ”Class segmentation and object localization with superpixel neighborhoods”, in _Proc. IEEE ICCV_, pp. 670–677, doi:[10.1109/ICCV.2009.5459175](http://dx.doi.org/10.1109/ICCV.2009.5459175).
# <a name=“GS09”></a>Grazzini J. and Soille P. (2009): [**Edge-preserving smoothing using a similarity measure in adaptive geodesic neighbourhoods**](http://www.sciencedirect.com/science/article/pii/S003132030800469X), _Pattern Recognition_, 42(10):2306-2316, doi:[10.1016/j.patcog.2008.11.004](http://dx.doi.org/10.1016/j.patcog.2008.11.004).
# <a name=“GDS10”></a>Grazzini J., Dillard S., and Soille P. (2010): [**Multichannel image regularisation using anisotropic geodesic filtering**](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5596008), in _Proc. IEEE ICPR_, pp. 2664-2667, doi:[10.1109/ICPR.2010.653](http://dx.doi.org/10.1109/ICPR.2010.653).
# <a name=“GS08”></a>Grazzini J.  and Soille P. (2008): [**Adaptive morphological filters using similarities based on geodesic time**](http://www.springerlink.com/content/f6v62233xqkklq72), in _Proc. DGCI_, Lecture Notes in Computer Science, vol. 4992, pp.519-528, doi:[10.1007/978-3-540-79126-3_46](http://dx.doi.org/10.1007/978-3-540-79126-3_46).
Hanbury, A., 2008. How do superpixels affect image segmentation? In: Ruiz-
Shulcloper, J., Kropatsch, W. (Eds.), Progress in Pattern Recognition, Image
Analysis and Applications - Proc. Iberoamerican Congress on Pattern Recognition.
No. 5197 in LNCS. Springer-Verlag, pp. 178–186.
# <a name=“LDM07”></a>R. Lerallut, E. Decenciere, and F. Meyer (2007): ”Image filtering using morphological amoebas”, _Image and Vision Computing_, 25(4):395–404, do:[10.1016/j.imavis.2006.04.018](http://dx.doi.org/10.1016/j.imavis.2006.04.018).
Levinshtein, A., Stere, A., Kutulakos, K., Fleet, D., Dickinson, S., Siddiqi, K.,
2009. TurboPixels: fast superpixels using geometric flows. IEEE Transactions
on Pattern Analysis and Machine Intelligence 31 (12), 2290–2297.
# <a name=“MPWMJ08”></a>A.P. Moore, S. Prince, J. Warrell, U. Mohammed, and G. Jones (2008): ”Superpixel lattices”, in _Proc. IEEE CVPR_, pp. 1–8, doi:[10.1109/CVPR.2008.4587471](http://dx.doi.org/10.1109/CVPR.2008.4587471).
# <a name=“FH04”></a>P. Felzenszwalb and D. Huttenlocher (2004): ”Efficient graph-based image segmentation”, _International Journal of Computer Vision_, 59(2):167-181, doi:[10.1023/B:VISI.0000022288.19776.77](http://dx.doi.org/10.1023/B:VISI.0000022288.19776.77).
Shi, J., Malik, J., 2000. Normalized cuts and image segmentation. IEEE Transactions
on Pattern Analysis and Machine Intelligence 22 (8), 888–905.
Soille, P., 2008. Constrained connectivity for hierarchical image decomposition
and simplification. IEEE Transactions on Pattern Analysis and Machine Intelligence
30 (7), 1132–1145.
