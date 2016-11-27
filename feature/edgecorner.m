%% EDGECORNER - Wrapper function for various edge/corner detection techniques.
%
%% Syntax
%     [Edge, Corner] = EDGECORNER(I);
%     [Edge, Corner] = EDGECORNER(I, met_edge, met_corner, 
%                                 'Property', propertyvalue, ...);
%
%% Inputs
% *|I|* : an input image with size |(X,Y,C), where |C>1| when |I| is multichannel.
%      
% *|met_edge|* : logical variable or string setting the method used for
%     detecting edges; it can be either:
%      
% * |'canny'| or |'log'| when performing Canny or log-Laplacian edge detection
%          [GW02] with function |EDGE|; in that case, the detector is applied
%          channel by channel,
% * |'vista'| for detecting edges likewise [PS06], a slightly variant of the
%          Canny detector implemented in Matlab,   
% * |'black'| for the anisotropic edge detector of [BSMH98] using the 
%          function |ANISOEDGE|, 
% * |'rothwell'| for the subpixel edge detector of [RMHN95] using the
%          function |ROTHWELLEDGE|, 
% * |'koethe'| for the tensor-based edge detector of [Koth03] using the
%          function |KOETHEDGE|, 
% * |'compass'| for the compass edge/corner detector of [RT01] using the
%          function |COMPASSEDGE|,
% * |'elder'| for the blur-based edge detector of [EZ98] using the function
%          |ELDERZUCKEREDGE|,
% * any of the strings |'prewitt'|, |'sobel'|, |'fast'|, |'derivative5'|,
%          |'derivative7'|, |'fleck'| or |'luengo'| for applying the 
%          corresponding edge detectors using the function |GRDSMOOTH|,
% * |'congrue'| for the edge/corner detection of [Koves03] based on
%          phase congruency and implemented in |CONGRUENCYEDGE|;
%      
% in the case |met_edge| is logical and |false|, then no edge are estimated,
%     it is expected that corners only are detected (see variable |met_corner| 
%     below); in the case |met_edge| is logical and |true|, it is assigned
%     its default value; default: |met_edge='vista'|.
%      
% *|met_corner|* : logical variable or string defining the method used for 
%     extracting corners; it can be either:
%      
% * |'harris'| or |'noble'| for Harris-Forster and Noble detections resp. 
%          using the function |HARRISCORNER|,
% * |'susan'| for Susan detection [SB97] with function |SUSANCORNER|,
% * |'fast9'|, |'fast10'|, |'fast11'|, |'fast12'|, for any FAST corner
%          detection using the function |FASTCORNER|,
% * |'cpda'| for curvature-based detection [ALFR09] using the function
%          |FASTCPDA|;
%      
% in the case |met_corner| is logical and |false|, then no corners are
%     extracted, it is expected that edges only are detected (see variable
%     |met_edge| above); in the case met_edge is logical and true, it is 
%     assigned its default value; default: |met_corner='harris'|.
%      
% Note that in some cases, the implemented edge detectors (see |met_edge|
% above) also enable to extract the corners (|'congrue'|, |'compass'| or
% |'koethe'|).
%
%% Property [propertyname  propertyvalues]
% *|'rho'|* : post-smoothing standard deviation, used either by either of the
%     detectors |'koethe'| and |'harris'| (see corresponding functions); 
%     default: |rho=1|.
%      
% *|'sig'|* : pre-smoothing standard deviation, used by various of the methods
%     above.
%      
% *|'reduce'|* : logical value or string defining the way the different channels
%     of a multispectral image are combined into the output edge map (see
%     |CANNYEDGE_BASE|); it can be either:
%      
% * |'igray'|: the input RGB image is converted to a gray image using 
%          function |RGB2GRAY|,
% * |'imax', 'isum'| : the input image (any dimension) is converted to a
%          gray image by taking the sum and the max over the different
%          channels resp.,
% * |'gmax'| : gradients are computed for the different channels and their
%          local pixelwise max is given as a single input to the edge
%          detector;
% * |'eor'| : calculations are made like for a multispectral image (as if
%          |reduce=false|), but the final edge map is taken as the logical 
%          |OR| of the output edge maps of the different channels;
%      
% in |true| case, it is set to |'sum'|; default: |reduce=false|, ie. no
%     'combination' is used, a multichannel map is output; this parameter
%     is naturally ignored when |I| is a scalar image; note that reduce
%     should be set to |true| when |met_edge='compass'|.
%      
% *|'int', 'samp'|* : parameters used for GST estimation in the method
%     |'koethe'|; see corresponding function and |GRD2GST|; default: 
%     |int='fast'| and |samp=1|.
% 
%% Outputs
% *|Edge|* : logical map of edges.
%      
% *|Corner|* : logical map of corners.
%      
% Note that if |met_edge| is set to |false| and only one lhs is required,
% then |Corner| is automatically passed in the output lhs variable (ie. output 
% variables' positions of |Edge| and |Corner| are exchanged).
%
%% References
% [Canny86]  J.F. Canny: "A computational approach to edge detection",
%      IEEE Trans. on Pattern Analysis and Machine Intelligence, 8(6):679-698,
%      1986.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4767851&tag=1>
%      
% [Zenzo86]  S. Di Zenzo: "A note on the gradient of a multi-image",
%      Computer Vision and Graphical Image Processing, 33:116-125, 1986.
%      <http://www.sciencedirect.com/science/article/pii/0734189X86902239>
%
% [HS88]  C.G. Harris and M.J. Stephens: "A combined corner and edge 
%      detector", Proc. Vision Conference, pp 147-151, 1988.
%      <http://www.bmva.org/bmvc/1988/avc-88-023.pdf>
%      
% [RMHN95] C. Rothwell, J. Mundy, B. Hoffman and V.-D. Nguyen: "Driving
%      Vision by Topology", Techn. report 2444, INRIA, 1995.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=477034&tag=1>
%      
% [HSSB97]  M. Heath, S. Sarkar, T. Sanocki, T. and K.W. Bowyer: "A robust
%      visual method for assessing the relative performance of edge-detection
%      algorithms, IEEE Trans. on Pattern Analysis and Machine Intelligence,
%      19(12):1338-1359, 1997.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=643893>
%      
% [SB97]  S.M. Smith and J.M. Brady: "{SUSAN} - A new approach to low  
%      level image processing", International Journal of Computer Vision,
%      23(1):45-78, 1997.
%      <http://www.lems.brown.edu/vision/courses/image-processing/Readings/smith95susan.pdf>
%      
% [EZ98]  J.H. Elder and S.W. Zucker: "Local scale control for edge 
%      detection and blur Estimation", IEEE Trans. on Pattern Analysis and
%      Machine Intelligence, 20(7), 1998.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=689301&tag=1>
%      
% [BSMH98] M. Black, G. Sapiro, D. Marimont and D. Heeger: "Robust 
%      anisotropic Diffusion", IEEE Trans. Image Processing, 7(3,):421-432, 
%      1998.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=661192>
%      
% [HSSB98]  M. Heath, S. Sarkar, T. Sanocki, and K. Bowyer: "Comparison of
%     edge detectors: a methodology and initial study",  Computer Vision and
%     Image Understanding, 69(1):38-54, 1998.
%     <http://www.sciencedirect.com/science/article/pii/S1077314297905877>
%
% [RT01]  M. Ruzon and C. Tomasi: "Edge, junction, and corner detection
%      using color distributions",  IEEE Transactions on Pattern Analysis
%      and Machine Intelligence, 23(11):1281-1295, 2001. 
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=969118>
%      
% [GW02]  R.C. Gonzales and R.E. Woods: "Digital Image Processing",      
%      Prentice Hall, 2002.
%      
% [Koves03]  P.D. Kovesi: "Phase congruency detects corners and edges",
%      Proc. DICTA, pp. 309-318, 2003. 
%      <www.csse.uwa.edu.au/~pk/research/pkpapers/phasecorners.pdf>
%      
% [Koth03]  U. Kothe: "Edge and junction detection with an improved 
%      structure tensor", Proc. of DAGM Symposium, LNCS 2781, pp. 25-32, 
%      Springer, 2003.
%      <http://hci.iwr.uni-heidelberg.de/Staff/ukoethe/papers/structureTensor.pdf>
%      
% [HY04]  X.C. He and N.H.C. Yung: "Curvature scale space corner detector
%      with adaptive threshold and dynamic region of support", Proc. ICPR,
%      vol. 2, pp. 791-794, 2004.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1334377>
%      
% [PS06]  L. Prasad and A. Skourikhine: "Vectorized image segmentation 
%      via trixel agglomeration", Pattern Recognition, 39:501-514, 2006.
%      <http://www.sciencedirect.com/science/article/pii/S0031320305003857>
%
% [ALFR09]  M. Awrangjeb, G. Lu, C.S. Fraser and M. Ravanbakhsh: "A fast
%      corner detector based on the chord-to-point distance accumulation 
%      Technique", Proc. DICTA, pp. 519-525, 2009.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5384897>
%
%% See also
% <EDGE.html |EDGE|>,
% <CORNER.html |CORNER|>,
% <CANNYEDGE.html |CANNYEDGE|>,
% <../../../vista/html/CANNYEDGES.html |CANNYEDGES|>,
% <CANNYEDGEPROD.html |CANNYEDGEPROD|>,
% <CANNYEDGEMAP.html |CANNYEDGEMAP|>,
% <SDGDEDGE.html |SDGDEDGE|>,
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>,
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <ANISOEDGE.html |ANISOEDGE|>,
% <ELDERZUCKEREDGE.html |ELDERZUCKEREDGE|>,
% <KOETHEDGE.html |KOETHEDGE|>,
% <PETROUEDGE.html |PETROUEDGE|>,
% <SUSANCORNER.html |SUSANCORNER|>,
% <HARRISCORNER.html |HARRISCORNER|>,
% <FASTCORNER.html |FASTCORNER|>,
% <FASTCPDA.html |FASTCPDA|>.
% Called:
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>.

%% Function implementation
function [Edge, Corner] = edgecorner(I, varargin)

%%
% parsing parameters

error(nargchk(1, 23, nargin, 'struct'));
error(nargoutchk(0, 2, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(I)
    error('edgecorner:inputerror','matrix required in input'); 
end

p = createParser('EDGECORNER');   
% principal optional parameters
p.addOptional('edge', 'vista', @(x) islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'canny','log','vista','koethe','black','rothwell', ...
    'sobel','prewitt','congrue','compass','elder', ...
    'fast','derivative5','derivative7','luengo','fleck'}))));
p.addOptional('corner', 'harris', @(x) islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'harris','noble','susan','cpda', ...
    'fast9','fast10','fast11','fast12'}))));
% other optional parameters
p.addParamValue('sig', 0.5, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addParamValue('rho', 1, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addParamValue('reduce', false, @(x)islogical(x) || ...
    (ischar(x) && any(strcmpi(x,{'igray','imax','isum','gmax','eor'}))));
% used by all corner detections
p.addParamValue('thres', 1, @(x)(isscalar(x) || numel(x)==2) && all(x>=0));
% used for method 'Koethe' only
p.addParamValue('int', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','conv','fast','ani'}))));
p.addParamValue('samp', 1, @(x)isscalar(x) && round(x)==x && x>=1 && x<=5);

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% check parameters

if islogical(p.edge) && ~p.edge && islogical(p.corner) && ~p.corner
    error('edgecorner:inputerror', ...
        'at least one of the variables met_corner or met_edge must be defined');
end

if islogical(p.edge) && p.edge
    p.edge = 'vista'; % reset to default
end    

if islogical(p.corner) && p.corner
    p.corner = 'harris';    
end

if strcmpi(p.edge,'compass') && size(I,3)==3 && ...
        ~(islogical(p.reduce) && p.reduce)
    warning('edgecorner:inputwarning', ...
        ['compass not using full color information like in original method - '...
        'set ''reduce'' to true for original approach']);
end

if strcmpi(p.edge,'congrue') && any(p.thres>1),
    warning('edgecorner:inputwarning', ...
        ['threshold should be in [0,1] ''with congrue'' method - '...
        'threshold set to 0.1']);
   p.thres = [0.1 0.1]; 
end

%% 
% main calculation

[Edge, Corner]  = edgecorner_base(I, p.edge, p.corner, ...
    p.rho, p.sig, p.thres, p.reduce, p.int, p.samp);

%%
% display

if p.disp
    nout = ~isempty(Edge) + ~isempty(Corner);
    figure,
    if ~isempty(Edge) 
        subplot(1,nout,1), imagesc(Edge), axis image off,
        title(['edge map - method ' p.edge]);
        if size(Edge,3)==1,  colormap gray;  end
    end
    if ~isempty(Corner)
        subplot(1,nout,1+~isempty(Edge)), imagesc(Corner), axis image off,
        title(['corner map - method ' p.corner]);
        if size(Corner,3)==1,  colormap gray;  end
    end
end

%%
% if only one lhs and Edge is empty, then exchange the variables
if nargout==1 && isempty(Edge)
    Edge = Corner;
    Corner = [];
end

end % end of edgecorner
