% grazzja - Utility function
%
% credit: J.Grazzini
% license: European Union Public License
%
% cluster/  Clustering and classification:		
%   fuzzbound      - Computes the fuzzy boundaries on a fuzzy categorical map.
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
%   bfs_base       - Compute breadth first search distances, times, and tree
%                    of a graph and extracts valid paths.
%   dfs_base       - Perform a depth-first search (DFS) of a graph.
%   scomponents_base - Compute the strongly connected components of a graph
%                    and extracts valid paths.
%   dijkstra       - Implementation of Dijkstra algorithm, allowing for single
%                    and multiple sources distance calculation. 
%   dijk_base      - Implementation of Dijkstra algorithm.
%   dijkadvanced_base - Calculate minimum costs and paths using Dijkstra's 
%                    algorithm.
%   graphmap       - Converts binary map to connected graph and reciprocally.
%   graph2map_base - Converts a spatial unweighted graph into a 2D logical
%                    map.
%   map2graph_base - Converts a 2D logical map into a weighted graph. 
%   graphweight_base - Convert an unweighted spatial graph into a weighted
%                    graph.
%   graphrepresent - Raster representation (and display) of a weighted graph. 
%   ixneighbours   - Return the indices of neighbour cells in a matrix. 
%
% pdem/  Pseudo-DEM generation functions:		
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
%   morphprofile   - Compute the morphological profile, the derivative 
%                    morphological profile by opening and closing and the
%                    morphological multiscale characteristics.
%   granulometry   - Compute a granulometry through a series of morphological
%                    operations.  
%   fractalmorph   - 
%   imreconstructby - Perform opening- an closing-by-reconstruction.
%   imrclose       - Wrapping function for closing by reconstruction.
%   imropen        - Wrapping function for opening by reconstruction.
%   imeragrad      - Compute the operation consisting of a morphological
%                    gradient followed by an erosion.
%   fastcrossdilation - 
%   imlabel        - Label connected components in 2-D arrays.  
%   bwisolated     - Extract isolated pixels from a logical image.
%   bwthinupsample - Perform the upsampled morphological thinning of a binary 
%                    image.
%   enhanceline    - 
%   flatstrel      - Wrapping function for STREL used for creating 2D flat
%                    structuring elements.
%
% propagation/  Propagation and Fast Marching methods:	
%   fmm            - Launch the Fast Marching algorithm in 2D.
%   fmmisopropagation_base - Apply classical and multistencil Fast Marching
%                    method in 2D.
%   fmmanisopropagation_base - Perform anisotropic Fast Marching.
%   fmmpath_base   - Extract a discrete geodesic path using gradient descent
%                    or Runge-Kutta method.
%   dijkstrapropagation_base - Shortest distance from multiple source points
%                    on graph.
%   im2front       - Perform a multiple front propagation from a set of seed
%                    points using a metric derived from the input image.
%   im2potential_base - Design a metric to be used as a potential function in
%                    front propagation.
%   potential2front_base - Propagate a (isotropic or anisotropic) front from
%                    a set of seed poinrs
%   saddlefront_base - Compute saddle points of an map of influence zones.
%
% segmentation/  Segmentation functions:	
%   geodesicwshed  - Perform geodesic based watershed segmentation of an image.
%   slicsuperpix   - Compute superpixels based on the Simple Linear Iterative
%                    Clustering segmentation technique.
%   geosuperpix    - Compute geodesic superpixels similarly to the Simple Linear
%                    Iterative Clustering segmentation technique.
%   lifetimeseg    - Pixelwise life time index computed over a stack of labelled
%                    images.
%   regionadjacency - Compute a sparse adjacency matrix for a labelled
%                    image of segmentation.
%   regionclean    - Clean up an image segmentation by deleting regions which
%                    do not obey certain shape/area criteria.
%
% statistics/  Basic statistical tools and operators:	
%   conv2fft       - FFT based two dimensional convolution.
%   convfft        - FFT based convolution and polynomial multiplication.
%   findzeroextrema1d - Wrapping function for popular Matlab functions 
%                    calculating extremas and/or zero-crossings of 1D signals. 
%   grididw_base   - Perform Inverse Distance Weighting or Simple Moving
%                    Average of sampled data.
%   grididwreduced_base - Perform reduced Inverse Distance Weighting.
%   histoequalization - Perform histogram equalization,.
%   integralimage_base - Compute the integral image of an image.
%   histointegral1d_base - Compute the integral histograms in cartesian space
%                    of an image.
%   histointegral2d_base - Compute the 2d integral histograms of the joint 
%                    distribution of a couple of images.
%   histointegraljoint_base - Compute the joint integral histograms.
%   rankintegral_base - 
%   histoslidefunction - Function for computing local statistical features
%                    based on local histograms computed over sliding windows.
%   localsum - Compute local sum of an image in a square window.
%   mutualsingfractal_base - Quantify the degree of mutual closeness of two
%                    singularity spectra.
%   accumarrayset  - Group elements from a data vector set as they share a
%                    same index; fast version of ACCUMARRAY(..., @(x){x}).
%   uniquenosort   - Find unique elements of vector likewise UNIQUE, but,
%                    contrary to it, without sorting the results.
%   allcombs       - Returns all 'crossed' combinations of elements in input
%                    vectors.
%   rescale        - Rescale data in [a,b].
%   clamp          - Clamp a value, cell of values or array.
%   nb_dims        - Debugged version of NDIMS.
%
% triangulation/  Triangulation segmentation and superpixel decomposition:
%   im2poly        - Decompose an image into significant polygons using the
%                    original VISTA segmentation based on constrained
%                    triangulation of the input image from edges. 
%   bound2vert     - Compute the so-called Attributed Contour Traces.
%   chainconnect_base - Produce contour chains of contiguous edges from a  
%                    edge map. 
%   chaincontour_base - Wrapping function for the original VISTA function
%                    contourchains.
%   chainlink_base - Links chains of edges using Kovesi's utility functions.
%   chaintrim_base - Trim the set of 2D points used for the constrained 
%                    Delaunay triangulation. 
%   vert2dtri      - Produce a constrained or unconstrained, Euclidean or non-
%                    Euclidean Delaunay triangulation from a set of vertices.
%   dtri2triedge   - 
%   dtriedgesample - Sampling of the underlying grid given a triangulation. 
%   dtriedgecircle_base - 
%   dtriedgecolor  - 
%   dtriedgecue - Simplifies the triangulation using spatial and spectral cues.
%   dtriedgecuespatial_base - 
%   dtriedgechord_base - 
%   dtriedgegestalt_base - 
%   dtriedgecuespectral_base - 
%   dtriedgecuetexture_base - 
%   dtriedge2poly - 
%   dtriedgepolydisplay - Function for displaying triangulations, trixels, 
%                    edges or polygons.
%
% texture/  Texture analysis and representations:		
%   localglcm2d    - Compute local textural features using the 2D histograms
%                    of joint distribution of greylevel pairs. 
%   localglov2d    - Compute local textural features using the 1D histograms
%                    of univariate distribution of greylevels.
%   localglsdv2d    - Compute local textural features using the 1D histograms
%                    of univariate distribution of greylevel differences.
%   integralglcm2d_base - 
%   integralglov2d_base - 
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
%   askifcontinue  - Utilitiy function to pause a program and asking if
%                    pursueing it or not.
%   catstruct      - Concatenate structures.
%   mergestruct    - Merge structures with unique fields.
%   str2vsubstr    - Vertically concatenate the substrings of a collection 
%                    of patterns strings found in the input string.
%   structa2astruct - Utility function converting a structure of arrays an 
%                    array of structures for utility.
%   comparecols    - Compare the columns of a matrix in increasing order.
%   createparser   - Create an instance of the inputParser class, with some
%                    additional fields.
%   getvarparser   - Returns the structure containing the full list of parameters
%                    from the input Parser structure.
%   freadenvi      - Read an image of ENVI standard type to a MATLAB array.

