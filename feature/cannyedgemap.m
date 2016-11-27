%% CANNYEDGEMAP - Find edge in an image from its directional derivatives.
%
%% Description
% This function finds edges in a image, given its directional derivatives,
% and possibly the magnitude/orientation images. It is nothing else than
% Canny's approach [Canny86] (non-maxima suppression + hystheresis
% thresholding) applied on already estimated derivatives.
%
%% Algorithm
% The Canny method finds edges by looking for local maxima of mag given as
% input. The method uses two thresholds, to detect strong and weak edges,
% and includes the weak edges in the output only if they are connected to 
% strong edges. This method is therefore less likely than the others to be 
% "fooled" by noise, and more likely to detect true weak edges.
%
%% Syntax
%     map = CANNYEDGEMAP(gx,gy);
%     map = CANNYEDGEMAP(gx,gy,'propertyname',propertyvalue,...)
%
%% Inputs
% *|gx, gy|* : directional derivatives of an image in X- (horizontal) and Y- 
%     (vertical) directions, typically estimated by differentiating the
%     Gaussian filter of an image; see |GRDSMOOTH| for discussion regarding
%     the direction/orientation of the gradient derivatives.
%
%% Property [propertyname  propertyvalues]
% *|'mag'|* : a matrix storing typically the image gradient magnitude or the
%     norm of the gradient structure tensor; if not parsed, it is the norm
%     of the |[gx,gy]| gradient field (and possibly taking the max over the
%     different channels if the gradient is 3D: see |CANNYEDGES|). 
%
% *|'or'|* : image containing feature normal orientation angles in degrees
%     (0-180), angles positive anti-clockwise.
%
% *|'der'|* : string setting the name of the function used for Canny edge
%     detection; it can be:
%
% * |'vista'| for the function |CANNYEDGES| used in Radius project,
% * |'edge'| (or |'matlab'|) for the standard |EDGE| detection
%          implemented in matlab (the same with |nbin=64| instead of 128 when 
%          estimating automatically the threshold values),
% * |'kovesi'| for the functions |NONMAXSUP| snd |HYSTTHRES| implemented 
%           by P.Kovesi [Kovesi];
%
% default: |der='edge'|.
%
% *|'hyst'|* : vector |[low high]| setting the low and high hystheresis threshold 
%     values with |0<=low<high<1| used with |der='edges'|; default: |hyst=[]|, 
%     and the threshold values are chosen automatically (see |EDGE|) or using
%     the variable |ratio| below.
%
% *|'ratio'|* : vector with pixel and threshold ratios used for automatically
%     selecting the hysteresis thresholds when those are not passed;
%     default: |ratio = [pixratio threshratio] = [1/3 0.08]|.
%
%% Output
% *|map|* : binary edge map with non maximally suppressed image.
%
%% Acknowledgment
% This function is a copy/paste and crop of several proposed functions 
% implemented for thresholding, non maximum suppression and thinning of
% edge maps.
%
%% References
% [Canny86]  J.F. Canny: "A computational approach to edge detection",
%      IEEE Trans. on Pattern Analysis and Machine Intelligence, 8(6):679-698,
%      1986.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4767851&tag=1>
% 
% [Kovesi]   P.D. Kovesi: "MATLAB and Octave Functions for Computer Vision
%      and Image Processing", The University of Western Australia.
%      <http://www.csse.uwa.edu.au/~pk/research/matlabfns/>
%
%% See also
% Related:
% <CANNYEDGE.html |CANNYEDGE|>,
% <EDGECORNER.html |EDGECORNER|>,
% <../../../vista/html/CANNYEDGES.html |CANNYEDGES|>,
% Called:
% <CANNYEDGEMAP_BASE.html |CANNYEDGEMAP_BASE|>. 

%% Function implementation
function map = cannyedgemap(gx,gy,varargin)

%%
% parsing and checking parameters

error(nargchk(1, 20, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

if ~(isnumeric(gx) && isnumeric(gy))
    error('cannyedgemap:inputparameter','matrices required in input'); 
elseif ~isequal(size(gx),size(gy))
    error('cannyedgemap:inputparameter', ...
        'input matrices with same size required in input'); 
end

p = createParser('CANNYEDGEMAP');   % create an instance of the inputParser class.
% function handling the function used for estimating Canny edges
p.addParamValue('mag', [],  @(x)isfloat(x));
p.addParamValue('or', [], @(x)isfloat(x));
p.addParamValue('der', 'matlab', @(x)ischar(x) && ...
    any(strcmpi(x,{'vista','edge','matlab','kovesi'})));
p.addParamValue('hyst', [], ...
    @(x)isvector(x) && length(x)<=2 && all(x>=0) && all(x<1));
p.addParamValue('ratio', [1/3 0.08], ...
    @(x)isvector(x) && length(x)==2 && all(x>=0) && all(x<1));

% parse and validate all input arguments
p.parse(varargin{:});
p = getvarParser(p);                                                            

%%
% main computation

map = cannyedgemap_base(gx, gy, p.der, p.mag, p.or, p.hyst, p.ratio);

%%
% display

if p.disp
    figure, imagesc(edgemap), axis image off, title('edge map');
    if size(edgemap,3) == 1, colormap 'gray'; end;
end

end  % end of cannyedgemap

