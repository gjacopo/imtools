**About**

Various [`Matlab`](https://nl.mathworks.com/) functions and tools used basic low-level image processing and analysis.

*version*:        1.0

*license*:      [EUPL](https://joinup.ec.europa.eu/sites/default/files/eupl1.1.-licence-en_0.pdf)

**Installation and usage**

Install the package `imtools` on a local drive and the paths of the different modules below to your `pathdef.m` setup file:
* `algebra`:  Matrix manipulation
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


See [`Contents.m`](Contents.m) for the whole list of functions. See [`Contents.html`](https://gjacopo.github.io/imtools/Contents.html) for a quick overview, 

That's it... Or almost: you may also want to install some useful `mex` built-in functions: run [`compile_feature.m`](feature/src/compile_feature.m),  [`compile_filter.m`](filter/src/compile_filter.m),  [`compile_propagation.m`](propagation/src/compile_propagation.m),  [`compile_pyramid.m`](pyramid/src/compile_pyramid.m) so as to compile some C/C++ source code. 

**Note** 
Run [`webpub`](https://gjacopo.github.io/imtools/misc/webpub.html) function to generate html documentation. Check script use.
