**About**

Various [`Matlab`](https://nl.mathworks.com/) functions and tools used basic low-level image processing and analysis.

*version*:        1.0

*license*:      [EUPL](https://joinup.ec.europa.eu/sites/default/files/eupl1.1.-licence-en_0.pdf)

**Installation and usage**

Install the package `imtools` on a local drive and the paths of the different modules below to your `pathdef.m` setup file:
* [`algebra`](algebra/index.md):  Matrix manipulation		
* [`derive`](derive/index.md):  Differentiation functions
* [`feature`](feature/index.md):  Feature extraction
* [`filter`](filter/index.md):  Filtering functions
* [`graph`](graph/index.md):  Graph and network analysis and representation
* [`surface`](surface/index.md):  Surface analysis functions
* [`pyramid`](pyramid/index.md):  Hierarchical pyramid analysis and decomposition	
* [`sharpen`](sharpen/index.md):  Sharpening and enhancing functions
* [`geometry`](geometry/index.md):  Basic geometrical tools
* [`kernel`](kernel/index.md):  Kernel operators
* [`morphology`](morphology/index.md):  Mathematical morphology
* [`propagation`](propagation/index.md):  Propagation and Fast Marching methods	
* [`segmentation`](segmentation/index.md):  Segmentation functions
* [`statistics`](statistics/index.md):  Basic statistical tools and operators
* [`texture`](texture/index.md):  Texture analysis and representations
* [`misc`](misc/index.md):  Miscellaneous functions* `algebra`:  Matrix manipulation

See [`Contents.m`](Contents.m) for the whole list of functions. See [`Contents.html`](https://gjacopo.github.io/imtools/Contents.html) for a quick overview, 

That's it... Or almost: you may also want to install some useful `mex` built-in functions: run [`compile_feature.m`](feature/src/compile_feature.m),  [`compile_filter.m`](filter/src/compile_filter.m),  [`compile_propagation.m`](propagation/src/compile_propagation.m),  [`compile_pyramid.m`](pyramid/src/compile_pyramid.m) so as to compile some C/C++ source code. 

**Note** 
Run [`webpub`](https://gjacopo.github.io/imtools/misc/webpub.html) function to generate html documentation. Check script use.
