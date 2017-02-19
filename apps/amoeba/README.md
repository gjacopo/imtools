amoeba
======

Amoeba-based superpixel partitioning of multispectral images into
---

**About**

Various basic low-level image processing tools used in the experiments/analyses reported in the references listed [below](References).

<table align="center">
    <tr> <td align="left"><i>documentation</i></td> <td align="left">available at: https://gjacopo.github.io/imtools/segmentation/amoebasuperpix.html</td> </tr> 
    <tr> <td align="left"><i>version</i></td> <td align="left">1.0 <i>(non-active development)</i> </td> </tr> 
    <tr> <td align="left"><i>since</i></td> <td align="left">2007</td> </tr> 
    <tr> <td align="left"><i>license</i></td> <td align="left"><a href="https://joinup.ec.europa.eu/sites/default/files/eupl1.1.-licence-en_0.pdfEUPL">EUPL</a>  <i>(cite the source code or any of the references below!)</i> </td> </tr> 
</table>

**Summary**

Segmentation is a difficult task because of the high complexity of images, where complexity refers to the large variety of pictorial representations of objects with the same semantic meaning and also to the extensive amount of available details.

We implement a new algorithm that:
# works essentially like a k-means based local clustering of pixels, but 
# enforces connectivity, so that it can efficiently generate compact, connected, and nearly uniform superpixels. 

This approach is based on the estimation of amoeba-like neighborhoods around selected cluster centers that exploit the connections between successive image pixels along geodesic paths in the image. The resulting superpixels capture the spatial/spectral redundancy in images and greatly reduce the complexity of subsequent image processing tasks. They provide convenient primitives from which to compute local image features when objects present in the scene have diverse scales or when they are not known in advance. 

<table>
<tr>
<td><kbd><img src="excerpt1.png" alt="doc SAS" width=“250"> </kbd></td>
<td><kbd><img src="exceprt1-amoeba-superpixels.png" alt="doc SAS" width=“250"> </kbd></td>
<td><kbd><img src="excerpt1-mean-amoeba-approximations.png" alt="doc R" width=“250"> </kbd></td>
</tr>
<header>
<td align="centre"><code>input excerpt</code></td>
<td align="centre"><code>amoeba superpixels (different scales displayed)</code></td>
<td align="centre"><code>amoeba mean-averaged approximations (different scales displayed)</code></td>
</header>
</table>

**Description**

The problem of segmenting an image into semantically meaningful units is a fundamental one in image processing. Indeed, visual information extraction is usually regarded as a segmentation issue [[FH04]](FH04). Besides practical issues (running time, user interaction, etc...), the key issue is what objects, if any, correspond to the segmented regions. It is well known that segmentation is an ill-posed problem whose ’correct’ solution is largely dependent on the application, if not completely subjective. Scene segmentation captures the local redundancy in the data, by reducing noise and variability, but aggregating pixels into segments often entails a decision that is unrelated to the final task. The goal is to perform this decision in a conservative way to minimize the risk of merging unrelated pixels from different objects. Over-segmentation is a pragmatic

**About**


**Description**

**Usage** 

Install Matlab package [`imtools`](https://gjacopo.github.io/imtools/) and run function  [`amoebasuperpix.m`](https://gjacopo.github.io/imtools/segmentation/amoebasuperpix.m).

**<a name="References"></a>References** 

# <a name=“ASSLFS10”></a>R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua and S. Susstrunk (2010): ”SLIC superpixels”, EPFL Technical Report no. 149300.
# <a name=“GDS10”></a>Grazzini J., Dillard S., and Soille P. (2010): [**Multichannel image regularisation using anisotropic geodesic filtering**](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5596008), in _Proc. IEEE ICPR_, pp. 2664-2667, doi:[10.1109/ICPR.2010.653](http://dx.doi.org/10.1109/ICPR.2010.653).
# <a name=“LDM07”></a>R. Lerallut, E. Decenciere, and F. Meyer (2007): ”Image filtering using morphological amoebas”, _Image and Vision Computing_, 25(4):395–404, do:[10.1016/j.imavis.2006.04.018](http://dx.doi.org/10.1016/j.imavis.2006.04.018).