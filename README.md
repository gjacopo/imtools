[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.322421.svg)](https://doi.org/10.5281/zenodo.322421)
imtools
=======

Tools (Matlab/C) for low-level image processing
---

**About**

Various basic low-level image processing tools used in the experiments/analyses reported in the references listed [below](#References).

<table align="center">
    <tr> <td align="left"><i>documentation</i></td> <td align="left">available at: https://gjacopo.github.io/imtools/</td> </tr> 
    <tr> <td align="left"><i>version</i></td> <td align="left">1.0 <i>(non-active development)</i> </td> </tr> 
    <tr> <td align="left"><i>since</i></td> <td align="left">2007</td> </tr> 
    <tr> <td align="left"><i>license</i></td> <td align="left"><a href="https://joinup.ec.europa.eu/sites/default/files/eupl1.1.-licence-en_0.pdfEUPL">EUPL</a>  <i>(cite the source code or any of the references below!)</i> </td> </tr> 
</table>

**Description**
Matlab tools for low-level image processing:
* `derive`:  Differentiation functions
* `feature`:  Feature extraction
* `filter`:  Filtering functions
* `graph`:  Graph and network analysis and representation
* `surface`:  Surface analysis functions
* `pyramid`:  Hierarchical pyramid analysis and decomposition		
* `sharpen`:  Sharpening and enhancing functions		 
* `geometry`:  Basic geometrical tools 	
* `kernel`:  Kernel operators
* `morphology`:  Mathematical morphology
* `propagation`:  Propagation and Fast Marching methods
* `segmentation`:  Segmentation functions	
* `statistics`:  Basic statistical tools and operators
* `texture`:  Texture analysis and representations		
* `misc`:  Miscellaneous functions	

See [`Contents.m`](Contents.m) for the whole list of functions. See file [`Contents.html`](https://gjacopo.github.io/imtools/Contents.html) for browsable version.

**Applications** 

Some examples of applications are described on the following pages:
* [amoeba superpixels](apps/amoebapix/) for superpixel partitioning of multispectral images into elementary, connected, nearly uniform units.
* [pseudo DEM](apps/geopdem/) for line networks extraction in digital images using geodesic propagation and flow simulation.

**<a name="References"></a>References** 

* Soille P. and Grazzini J. (2011): [**IMAGE-2006 Mosaic: Geometric and radiometric consistency of input imagery**](http://publications.jrc.ec.europa.eu/repository/bitstream/JRC49168/lbne23636enn.pdf), _Publications Office of the European Union_, doi:[10.2788/50967](http://dx.doi.org/10.2788/50967).
* Grazzini J., Dillard S., and Soille P. (2010): [**A new generic method for the semi-automatic extraction of river and road networks in low and mid-resolution satellite images**](http://spiedigitallibrary.org/proceedings/resource/2/psisdg/7830/1/783007_1), in _Proc. of SPIE_ - Image and Signal Processing for Remote Sensing XVI, vol. 7830, pp. 7830071-10, doi:[10.1117/12.865052](http://dx.doi.org/10.1117/12.865052).
* Grazzini J., Dillard S., and Prasad L. (2010): [**Simultaneous hierarchical segmentation and vectorization of satellite images through combined non-uniform data sampling and anisotropic triangulation**](https://spie.org/Publications/Proceedings/Paper/10.1117/12.865047), in _Proc. of SPIE_ - Image and Signal Processing for Remote Sensing XVI, vol. 7830, pp. 78300F-13, doi:[10.1117/12.865047](http://dx.doi.org/10.1117/12.865047).
* Dillard S., Prasad L., and Grazzini J. (2010): [**Region and edge-adaptive sampling and boundary completion for segmentation**](http://link.springer.com/chapter/10.1007/978-3-642-17274-8_7), in _Proc AVC_, Lecture Notes in Computer Science, vol. 6454, pp.64-74, doi:[10.1007/978-3-642-17274-8_7](http://dx.doi.org/10.1007/978-3-642-17274-8_7).
* Grazzini J., Dillard S., and Soille P. (2010): [**Multichannel image regularisation using anisotropic geodesic filtering**](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5596008), in _Proc. IEEE ICPR_, pp. 2664-2667, doi:[10.1109/ICPR.2010.653](http://dx.doi.org/10.1109/ICPR.2010.653).
* Grazzini J. and Soille P. (2010): [**Iterative ramp sharpening for structure/signature-preserving simplification of images**](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5597348), in _Proc. IEEE ICPR_, pp. 4586-4589, doi:[10.1109/ICPR.2010.1120](http://dx.doi.org/10.1109/ICPR.2010.1120).
* Grazzini J.  and Soille P. (2009): [**Image filtering based on locally estimated geodesic functions**](http://www.springerlink.com/content/v264v11754004500), in _Proc. VISIGRAPP_, Communications in Computer and Information Science, vol. 24, pp.123-124, doi:[10.1007/978-3-642-10226-4_10](http://dx.doi.org/10.1007/978-3-642-10226-4_10).
* Soille P. and Grazzini J. (2009): [**Constrained connectivity and transition regions**](http://www.springerlink.com/content/g6h8mk8447041532/), in _Proc. ISMM_, Lecture Notes in Computer Science, vol. 5720, pp.59-69, doi:[10.1007/978-3-642-03613-2_6](http://dx.doi.org/10.1007/978-3-642-03613-2_6).
* Grazzini J. and Soille P. (2009): [**Edge-preserving smoothing using a similarity measure in adaptive geodesic neighbourhoods**](http://www.sciencedirect.com/science/article/pii/S003132030800469X), _Pattern Recognition_, 42(10):2306-2316, doi:[10.1016/j.patcog.2008.11.004](http://dx.doi.org/10.1016/j.patcog.2008.11.004).
* Grazzini J.  and Soille P. (2008): [**Adaptive morphological filters using similarities based on geodesic time**](http://www.springerlink.com/content/f6v62233xqkklq72), in _Proc. DGCI_, Lecture Notes in Computer Science, vol. 4992, pp.519-528, doi:[10.1007/978-3-540-79126-3_46](http://dx.doi.org/10.1007/978-3-540-79126-3_46).
* Grazzini J., Soille P., and Bielski C. (2007): [**On the use of geodesic distances for spatial interpolation**](http://www.geocomputation.org/2007/7C-Spatial_statistics_3/7C3.pdf), in _Proc GeoComputation_, 7C3, uri:http://publications.jrc.ec.europa.eu/repository/handle/JRC37427.
* Bielski C., Grazzini J., and Soille P. (2007): [**The Little algorithm that grew: Scaling the morphological image compositing algorithm to meet the chanllenges of processing large image data sets**](http://www.geocomputation.org/2007/1A-Remote_Sensing_1/1A5.pdf), in _Proc GeoComputation_, 7C3, uri:http://publications.jrc.ec.europa.eu/repository/handle/JRC44083.
* Bielski C., Grazzini J., and Soille P. (2007): [**Automated morphological image composition for mosaicing large image data sets**](http://ieeexplore.ieee.org/document/4423743/), in _Proc IEEE IGARSS_, doi:[10.1109/IGARSS.2007.4423743](http://dx.doi.org/10.1109/IGARSS.2007.4423743).
* Grazzini J. and Soille P. (2007): [**Improved morphological interpolation of elevation contour data with generalised geodesic propagations**](http://link.springer.com/chapter/10.1007%2F978-3-540-74272-2_92), in _Proc. CAIP_, Lecture Notes in Computer Science, vol. 4673, pp.742-750, doi:[10.1007/978-3-540-74272-2_92](http://dx.doi.org/10.1007/978-3-540-74272-2_92).
* Soille P. and Grazzini J. (2007): [**Extraction of river networks from satellite images by combining morphology and hydrology**](http://link.springer.com/chapter/10.1007%2F978-3-540-74272-2_79), in _Proc. CAIP_, Lecture Notes in Computer Science, vol. 4673, pp.636-644, doi:[10.1007/978-3-540-74272-2_79](http://dx.doi.org/10.1007/978-3-540-74272-2_79).
