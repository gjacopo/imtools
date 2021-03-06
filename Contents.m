% grazzja - Utility function
%
% credit: J.Grazzini
% license: European Union Public License
%
% algebra/  Matrix manipulation:
%   accumarrayset  - Overload ACCUMARRAY to group elements from a data vector set.
%   allcombs       - Retrieve all 'crossed' combinations of elements in various 
%                    input vectors.
%   allvcombs      - Retrieve all 'crossed' combinations of elements of a single
%                    input vector.
%   cellnumsubtrim - Trim a set of numerical sequences from its subsequences.
%   clamp          - Clamp data.
%   compcol        - Compare columns of a matrix.
%   ind2rc         - Variant of IND2SUB behaving nicely with single output.
%   mat2rc         - Convert a matrix to a row or column vector using the COLON
%                    operator.
%   nb_dims        - Overload NDIMS.
%   posinsert      - Insert numeric values in a vector or matrix at given positions.
%   rescale        - Rescale data.
%   reverse        - Flip a vector.
%   triunpack      - Unpack the non-redundant representation of a triangular or 
%                    symmetric matrix..
%   uniqueunsort   - Retrieve unsorted list of unique elements.  
%
% derive/  Differentiation functions:		
%   gstsmooth      - Compute the Gradient Structure Tensor of a (possibly 
%                    multichannel) image.
%   grdsmooth      - Compute the directional derivatives of an image.
%   hessmooth      - Compute the Hessian Tensor of a (possibly multichannel)
%                    image.
%   gstdecomp      - Perform the eigen decomposition or the reconstruction
%                    of a tensor field. 
%   gstfeature     - Derive features like norm, orientation, coherence, or 
%                    inertia from tensor representation
%   grd2gst        - Compute the entries of the gradient structure tensor of
%                    an image from its directional derivatives.
%   grd2hess       - Compute the entries of the Hessian of an image from its
%                    directional derivatives.
%   grdmask        - Compute the directional derivatives of an image through
%                    classical gradient operators using local difference masks.
%
% feature/  Feature extraction:		
%   anisoedge      - Anisotropic diffusion edge detector.
%   cannyedge      - Wrapping function calling either edge, cannyedges or
%                    cannyedgemap.
%   cannyedgefuse  - Two-scale Canny edge detection by fusion of the binary
%                    Canny edge maps at two different scales.
%   cannyedgemap   - Edge map calculation. 
%   cannyedgeprod  - Two-scale Canny edge detection obtained by multiplication
%                    of Canny edge responses at two different scales.
%   compassedge    - Wrapping function for the Generalized Compass Operator.
%   congruencyedge - Wrapping function for Kovesi'sedge and corner detector
%                    based on phase congrunecy. 
%   corner         - Wrapping function for various corner detection methods.
%   curvecorner    - Corner detection in edge map based on curvature. 
%   edgecorner     - Wrapping function for various edge/corner detection
%                    techniques.
%   elderzuckeredge - Elder & Zucker's multi-scale edge detector.
%   fastcorner     - Wrapping function for the FAST corner detector.
%   fastcpda       - Fast corner detector based on the chord-to-point distance
%                    accumulation technique.
%   fuzzbound      - Compute the fuzzy boundaries on a fuzzy categorical map.
%   harriscorner   - Harris' corner/keypoints detector.
%   inforgbedge    - Provide edge magnitude and orientation of colour images
%                    based on its derivatives
%   koethedge      - Koethe's edge/corner detector based on the gradient
%                    structure tensor.
%   multiedge      - Compute a multiscale edge map using a 2D Sobel-like
%                    wavelet.
%   petrouedge     - Petrou & Kitller's edge detector.
%   rothwelledge   - Rothwell's subpixel edge detector.
%   sdgdedge       - Edge detection filter based on the 2nd derivative and/or
%                    Laplacian in the gradient direction.  
%   susancorner    - SUSAN Edge/corner detection or image smoothing.
%
% filter/  Filtering functions:
%   adaptivefilt   - Perform adaptive filtering.
%   anisor         - Compute orientation information from the anisotropic
%                    continous scale model.
%   anisorfilt     - Morphological anisotropic diffusion using the orientation
%                    information from the anisotropic continous scale model.
%   bilateral_base - Perform bilateral filtering on monochrome images. 
%   bilateralstack_base - Fast version of the bilateral filter using stacks.
%   blurmap        - Adaptive local isotropic Gaussian bluring of an image.
%   convolution    - Compute convolution with centered filter.
%   findlocalextrema - Extracts the extrema of any multispectral image in
%                    local square neigbourhoods.
%   findlocalmax   - Determine the local maxima of a given matrix image.
%   fuzzimpulse    - Perform fuzzy impulse noise reduction in images.
%   geodesicfilt   - Sigma-filter based on local geodesic similarity measures 
%                    for edge-preserving smoothing of multispectral images.
%   mdlfilt        - Local adaptive filtering based on Minimal Description
%                    Length principle.
%   smoothfilt     - Perform either classical isotropic (Gaussian) smoothing
%                    or anisotropic non-linear  smoothing.
%   tensanifilt    - Perform Gaussian adaptive filtering of an image using
%                    an anisotropic tensor field.
%   tenscaledirfilt - Perform Gaussian adaptive directional filtering along
%                    tensor field using a local adaptive scale.
%   wmedfilt2      - Perform weighted median filtering of a (possibly
%                    multispectral) image.
%   wordfilt2      - Perform weighted order-statistic filtering of a
%                    (possibly multispectral) image.
% 
% fractal/  Fractal analysis:
%   fractalwave     - 
%   fractalwavestat - 
%   singfractalrecons_base - 
%   waveprofile_base - 
%
% graph/  Graph and network analysis and representation:
%   bfs            - Compute breadth first search distances, times, and tree
%                    of a graph and extracts valid paths.
%   dfs            - Perform a depth-first search (DFS) of a graph.
%   dijk           - Implementation of Dijkstra algorithm.
%   dijkstra       - Implementation of Dijkstra algorithm, allowing for single
%                    and multiple sources distance calculation. 
%   dijkadvanced_base - Calculate minimum costs and paths using Dijkstra's 
%                    algorithm.
%   graphmap       - Convert binary map to connected graph and reciprocally.
%   pathfinder     - Find paths in a graph.
%   graphweight    - Convert an unweighted spatial graph into a weighted graph.
%   graph2map_base - Convert a spatial unweighted graph into a 2D logical
%                    map.
%   map2graph_base - Convert a 2D logical map into a weighted graph. 
%   graphweight_base - Convert an unweighted spatial graph into a weighted
%                    graph.
%   graphdisplay   - Raster representation (and display) of a weighted graph.
%   ixneighbours   - Return the indices of neighbour cells in a matrix. 
%   scomponents    - Compute the strongly connected components of a graph
%                    and extracts valid paths.
%   triadjacency   - Construct the adjacency matrix associated to an order-3 graph.
%
% surface/  Surface analysis functions:		
%   pdem           - Pseudo DEM generation from a singular optical image.
%   wflowacc       - Upslope area (flow accumulation) algorithm for Digital 
%                    Elevation Models.
%   sflowinf       - Upslope area algorithm for Digital Elevation Models using 
%                    single flow direction dinf.
%   sflow8         - Upslope area algorithm for Digital Elevation Models using 
%                    single flow direction d8 and gd8.
%   sfdacc  - Upslope area algorithm for Digital Elevation Models using single 
%                    flow direction.
%
% pyramid/  Hierarchical pyramid analysis and decomposition:		
%   pyrdown        - Top-down pyramidal hierarchical decomposition.
%   pyrup          - Bottom-up pyramidal hierarchical decomposition.  
%   reduce2d       - Wrapping function for the reduce2D operator of the BIG
%                    Spline Pyramids Software.
%   expand2d       - Wrapping function for the expand2d operator of the BIG
%                    Spline Pyramids Software.
%
% sharpen/  Sharpening and enhancing functions:		 
%   rampsharp       - Iterative ramp sharpening of multispectral images.
%   maptransition   - Transition pixels mapping.
%   localorientzone - Define the clockwise zone of an image of orientation.
%   localorientfeature - Local zone-based filtering.
%
% geometry/  Basic geometrical tools: 	
%   bresenhamline  - 2D Bresenham line algorithm. 
%   downscalexy    - Downsample monotically an input image. 
%   upscalexy      - Upsample monotically an input image.
%   gridblk        - Compute the centers' coordinates and the size of the
%                    square blocks dividing an image in regular tiles.
%   sampletriangle - Sample lattice points inside triangles.
%   sampledisk     - Sample lattice points inside disks.
%   drawleakcircle - Create and/or draw a discrete disk with holes. 
%   drawcircle     - Create and/or draw a simple circle.
%   drawgestalt    - Create and/or draw simple figures of Gestalt contours.
%
% kernel/  Kernel operators:		
%   gausskernel    - Compute a 1D or 2D Gaussian filter.
%   dirgausskernel - Compute an oriented Gaussian kernel.
%   hourglasskernel - Compute an oriented non-linear spatial hourglass filter.
%   euclidkernel   - Create a kernel spatially weighted with the euclidean
%                    distance.
%   local3x3kernel - Compute local kernels defined by zones and levels.
%   neiposkernel   - Create the kernel containing the indices and subscript
%                    positions of the grid pixels of a given square window.
%
% morphology/  Mathematical morphology:	
%   asf            - Apply the Alternate Sequential Filtering and the ASF by 
%                    reconstruction of an image.
%   morphprofile   - Compute morphological profiles by opening and closing and 
%                    the morphological multiscale characteristics.
%   granulometry   - Compute a granulometry through a series of morphological 
%                    operations.
%   fractalmorph   - Fractal morphological decomposition.
%   imreconstructby - Perform opening- an closing-by-reconstruction.
%   imrclose       - Wrapping function for closing by reconstruction.
%   imropen        - Wrapping function for opening by reconstruction.
%   imeragrad      - Compute the operation consisting of a morphological gradient 
%                    followed by an erosion.
%   imlabel        - Label connected components in 2-D arrays.
%   bwisolated     - Extract isolated pixels from a logical image.
%   bwshrink1D     - Shrinking of 1D binary signals.
%   bwthinupsample - Perform the upsampled morphological thinning of a binary image.
%   flatstrel      - Wrapping function for STREL used for creating 2D flat structuring 
%                    elements.
%
% propagation/  Propagation and Fast Marching methods:	
%   fmm            - Launch the Fast Marching algorithm in 2D.
%   fmmisopropagation - Apply classical and multistencil Fast Marching method in 2D.
%   fmmanisopropagation - Perform anisotropic Fast Marching.
%   fmmpath   - Extract a discrete geodesic path using gradient descent or 
%                    Runge-Kutta method.
%   dijkstrapropagation - Shortest distance from multiple source points on graph.
%   im2front       - Perform a multiple front propagation from a seed points 
%                    using a metric derived from the input image.
%   im2potential - Design a metric to be used as a potential function in front 
%                    propagation.
%   potential2front - Propagate a (isotropic or anisotropic) front from a set 
%                    of seed points.
%   saddlefront_base - Compute saddle points of an map of influence zones.
%   im2amoeba - Amoeba-like front propagation over an image.
%
% segmentation/  Segmentation functions:	
%   amoebasuperpix - Amoeba-like superpixel segmentation.
%   slicsuperpix   - Compute superpixels based on the Simple Linear Iterative
%                    Clustering segmentation technique.
%   geosuperpix    - Compute geodesic superpixels similarly to the Simple Linear
%                    Iterative Clustering segmentation technique.
%   geodesicwshed  - Perform geodesic based watershed segmentation of an image.
%   lifetimeseg    - Pixelwise life time index computed over a stack of labelled
%                    images.
%   regionadjacency - Compute a sparse adjacency matrix for a labelled
%                    image of segmentation.
%   regionclean    - Clean up an image segmentation by deleting regions which
%                    do not obey certain shape/area criteria.
%
% statistics/  Basic statistical tools and operators:	
%   compressrange           - Compress the entries of an image.
%   findimportantextrema1D - Find important extrema of 1D signals.
%   findzeroextrema1d - Wrapping function for popular Matlab functions 
%                    calculating extremas and/or zero-crossings of 1D signals. 
%   grididw_base   - Perform Inverse Distance Weighting or Simple Moving
%                    Average of sampled data.
%   grididwreduced_base - Perform reduced Inverse Distance Weighting.
%   histoequalization - Perform histogram equalization,.
%   integralimage - Compute the integral image of an image.
%   integralhisto1d - Compute the integral histograms in cartesian space of an image.
%   integralhisto2d - Compute the 2d integral histograms of the joint distribution of a 
%                    couple of images.
%   integralhistojoint - Compute the joint integral histograms.
%   integralrank - Compute rank image.

%   slidehistofun - Function for computing local statistical features
%                    based on local histograms computed over sliding windows.
%   localsum - Compute local sum of an image in a square window.
%   mutualsingfractal_base - Quantify the degree of mutual closeness of two
%                    singularity spectra.
%   mutualspectra - Measure the mutual closeness of singularity spectra.
%   stats3x3 - Basic 2D statistical operations in local 3x3 neighbourhoods.
%   pqueue - Implement a priority queues.
%   nanaverage - Overload NANMEAN.
%   nansummation - Overload NANSUM.
%
% texture/  Texture analysis and representations:		
%   localglcm2d    - Compute local textural features using the 2D histograms
%                    of joint distribution of greylevel pairs. 
%   localglov2d    - Compute local textural features using the 1D histograms
%                    of univariate distribution of greylevels.
%   localglsdv2d    - Compute local textural features using the 1D histograms
%                    of univariate distribution of greylevel differences.
%   histfeatures   - Calculate all textural features from (co)occurrences
%                    matrices.
%   histcontrast   - Calculate the contrast feature .
%   histcorrelation - Calculate the correlation feature.
%   histdissimilarity - Calculate the dissimilarity feature.
%   histenergy     - Calculate the ENERGY (2nd angular moment) feature.
%   histentropy10  - Calculate the entropy feature.
%   histentropy2   - Calculate the entropy feature in bits.
%   histhomogeneity - Calculate the homogeneity feature.
%   histidifference - Calculate the inverse difference feature.
%   histdissimilarity - Calculate the maximum feature.
%   histmean       - Calculate the mean feature.
%   histvariance   - Calculate the variance feature.
%
% misc/  Miscellaneous functions:		
%   askifcontinue  - Utilitiy function to pause a program execution.
%   catstruct      - Concatenate structures.
%   cellstrsubtrim - Trim a set of strings from its matching substrings.
%   createparser   - Create an instance of the inputParser class, with some 
%                    additional fields.
%   editemplate    - Create and edit a new function with standard documentation 
%                    template.
%   getvarparser   - Return the structure containing the full list of parameters 
%                    from the input Parser structure.
%   legendgroup    - Provide legends for groups of plot handles.
%   mergestruct    - Merge structures with unique fields.
%   readenviraster - Read an image of ENVI standard type to a MATLAB array.
%   readenviroi    - Read ENVI files.
%   str2vsubstr    - Vertically concatenate the substrings of a collection of 
%                    patterns strings found in the input string.
%   structa2astruct - Convert a structure of arrays to an array of structures.
%   suptitle       - Utility to put a title above all subplots
%   webpub         - Creation of web formatted functions' documentation.
%   whichpath      - Return the full path name of a function.


