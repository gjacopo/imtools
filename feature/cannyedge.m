%% CANNYEDGE - Wrapping function for EDGE, CANNYEDGES and CANNYEDGEMAP.
% 
%% Description
% Nothing else than a wrapping function calling either |EDGE|, |CANNYEDGES| 
% or  |CANNYEDGEMAP|, written for convenience (and consistency with the
% other edge detection methods implemented in functions |EDGE<METHOD>|).
%
%% Syntax
%     edgemap = CANNYEDGE(I,sigma);
%     [edgemap, mag, or] = CANNYEDGE(I, sigma, ...
%                               'Property', propertyvalue, ...);
%
%% Inputs
% *|I|* : input image.
%
% *|sigma|* : optional standard deviation of the Gaussian filter used for
%     smoothing of the image; default: |sigma=1|.
%
%% Property [propertyname  propertyvalues]
% *|'der'|* : optional string storing the method used: |'edge'| (or |'matlab'|),
%     |'vista'|, |'kovesi'|, or any of the other options used by |GRDSMOOTH|:
%     |'fast'|, |'conv'|, |'fleck'|, |'opt'|, |'tap5'|, |'tap7'|, |'sob'|, 
%     |'prew'|, |'opt'|, |'circ'|, |'ana'| or |'lue'|; default: |der='edge'|.
%
% *|'reduce'|* : logical value or string defining the way the different channels
%     of a multispectral image are combined into the output edge map (see
%     |CANNYEDGE_BASE|); it can be either:
%
% * |'igray'| : the input RGB image is converted to a gray image using the
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
%     is naturally ignored when |I| is a scalar image.
%
% *|'hyst'|* : vector |[low high]| setting low and high hystheresis threshold 
%     values with |0<=low<high<1| used with |der='edge'|; default: |hyst=[]|, 
%     and the threshold values are chosen automatically (see |EDGE|).
%
% *|'samp'|* : to perform |(* samp)| interpolation of the input image to avoid
%     aliasing; default: |samp=1|.
%
%% Outputs
% *|edgemap|* : logical map of detected edges.
%
% *|mag, or|* : optional output matrices storing the magnitude and orientation
%     of the estimated gradient field.
%
%% Reference
% [Canny86]  J.F. Canny: "A computational approach to edge detection",
%      IEEE Trans. on Pattern Analysis and Machine Intelligence, 8(6):679-698,
%      1986.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4767851&tag=1>
%
%% See also
% Related:
% <CANNYEDGES.html |CANNYEDGES|>,
% <EDGECORNER.html |EDGECORNER|>,
% <CANNYEDGEPROD.html |CANNYEDGEPROD|>,
% <ROTHWELLEDGE.html |ROTHWELLEDGE|>,
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>,
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <ANISOEDGE.html |ANISOEDGE|>,
% <ELDERZUCKEREDGE.html |ELDERZUCKEREDGE|>,
% <KOETHEDGE.html |KOETHEDGE|>,
% <SDGDEDGE.html |SDGDEDGE|>,
% <PETROUEDGE.html |PETROUEDGE|>,
% <CANNYEDGEPROD.html |CANNYEDGEPROD|>,
% <../../derive/html/GRDSMOOTH_BASE.html |GRDSMOOTH_BASE|>.
% Called:
% <CANNYEDGE_BASE.html |CANNYEDGE_BASE|>. 

%% Function implementation
function [edgemap,varargout] = cannyedge(I,varargin)

%% 
% parsing parameters

error(nargchk(1, 18, nargin, 'struct'));
error(nargoutchk(1, 3, nargout, 'struct'));

if ~isnumeric(I)
    error('cannyedge:inputparameter','a matrix is required in input'); 
end

p = createParser('CANNYEDGE');   
p.addOptional('sigma', 1, @(x)isscalar(x) && x>=0.05);
p.addParamValue('der', 'matlab', @(x)ischar(x) && ...
    any(strcmpi(x,{'edge','matlab','vista','kovesi','fast','conv','fleck', ...
    'opt','tap5','tap7','sob','prew','opt','circ','ana','lue'})));
p.addParamValue('reduce', false, @(x)islogical(x) || ...
    (ischar(x) && any(strcmpi(x,{'igray','imax','isum','gmax','eor'}))));
p.addParamValue('hyst', [], @(x)isempty(x) || ...
    (@(x)isvector(x) && length(x)<=2 && all(x>=0) && all(x<1)));
p.addParamValue('samp', 1, @(x)isscalar(x) && round(x)==x && x>=1 && x<=5);
% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% checking/setting variables
C = size(I,3);

if C==1 && ~(islogical(p.reduce) && ~p.reduce), 
    warning('cannyedge:inputwarning', ...
        'ignored option ''reduce'' with scalar image');
    p.reduce = false;  
end

if C~=3 && ischar(p.reduce) && strcmpi(p.reduce,'igray')
    error('cannyedge:inputerror', ...
         ['option ''reduce'' set to ' p.reduce ' for RGB images only']);
end

if strcmpi(p.der,'matlab') 
    if nargout>=2 
    warning('cannyedge:inputwarning', ...
        'no magnitude/orientation computed with method ''matlab''');
    varargout = [];
    elseif ischar(p.reduce) && strcmpi(p.reduce,'gmax')
        warning('cannyedge:inputwarning', ...
            'option ''gmax'' ignored with method ''matlab''');
        p.reduce = false;
    end
end

%%
% main calculation

if nargout==1
    edgemap = cannyedge_base(I, p.sigma, p.der, p.samp, p.hyst, p.reduce);

else
    [edgemap, mag, or] = cannyedge_base(I, p.sigma, p.der, p.samp, p.hyst, p.reduce);
    if nargout>=2,     varargout{1} = mag;     end
    if nargout==3,     varargout{2} = or;     end
    
end

%%
% display

if p.disp
    figure, imagesc(edgemap), axis image off, title('Canny edge map');
    if size(edgemap,3) == 1, colormap 'gray'; end;
    if nargout==3
        figure, subplot(1,2,1)
        imshow(mag), axis image off, title('gradient magnitude');
        subplot(1,2,2)
        imshow(rescale(or,0,1)), axis image off, 
        title('gradient orientation');
        
    end
end

end % end of cannyedge
