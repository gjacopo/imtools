%% CANNYEDGEFUSE - Two-scale Canny edge detection.
%
%% Description
% Perform two-scale Canny edge detection through (logical) multiplication of 
% the binary edge maps calculated using two Canny edge detectors at different
% scales.
%
%% Syntax
%     edge = CANNYEDGEFUSE(I, sig2);
%     edge = CANNYEDGEFUSE(I, sig2, sig1, ser, ...
%                      'propertyname',propertyvalue,... );
% 
%% Inputs
% *|I|* : input image.
%      
% *|sig2, sig1|* : standard-deviations of, resp., the coarsest and the 
%     finest Gaussian filters used for edge detection; default: |sig2| needs
%     always to be passed as input and |sig1=1|.
%      
% *|serad|* : radius of the disk-shaped structuring element used for dilating
%     the edge map obtained at the coarsest scale; default: |serad=3|.
%
%% Property [propertyname  propertyvalues]
% *|'der'|* : string setting the name of the function used for Canny edge
%     detection; it can be:
%      
% * |'vista'| for the function |CANNYEDGES| used in Radius project,
% * |'matlab'| or |'edge'| for the standard |EDGE| detection implemented in
%          matlab (the same with |nbin=64| instead of 128 when estimating
%          automatically the threshold values),
% * |'kovesi'| for the function |CANNY| implemented by P.Kovesi;
%      
% default: |der='matlab'|; see |CANNYEDGEMAP|.
%      
% *|'hyst'|* : vector |[low high]| setting low and high hystheresis threshold 
%       values with |0<=low<high<1| used with |der='edge'|; default: |hyst=[]|, 
%       and the threshold values are chosen automatically (see |EDGE|).
%      
% *|'reduce'|* : logical value or string defining the way the different channels
%       of a multispectral image into the edge map (see |CANNYEDGE_BASE|); it 
%       can be either:
%      
% * |'igray'|: the input RGB image |I| is converted to a gray image using the
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
% in true case, it is set to |'sum'|; default: |reduce=false|, ie. no
%     'combination' is used, a multichannel map is output; this parameter
%     is naturally ignored when |I| is a scalar image.
%
%% Output
% *|edgemap|* : a binary edge map.
%
%% Reference
% [PCPN06]  G. Papari P. Campisi, N. Petkov and A. Neri: "A multiscale 
%      approach to contour detection by texture suppression", Proc. SPIE 
%      2006 Image processing: Algorithm and Systems, vol. 6064A, 2006.
%      <http://spiedigitallibrary.org/proceedings/resource/2/psisdg/6064/1/60640D_1>
% 
%% See also
% Called: 
% <CANNYEDGEFUSE_BASE.html |CANNYEDGEFUSE_BASE|>.

%% Function implementation
function edgemap = cannyedgefuse(I, sig2, varargin)

%% 
% parsing and checking parameters

error(nargchk(1, 18, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

if ~isnumeric(I)
    error('cannyedgefuse:inputparameter','a matrix is required in input'); 
end

p = createParser('CANNYEDGEFUSE'); 
% function handling the function used for estimating Canny edges
p.addRequired('sig2', @(x)x>0.15);
p.addOptional('sig1', 1., @(x)isscalar(x) && isfloat(x) && x>0);
p.addOptional('serad', 3, @(x)isscalar(x) && isinteger(x) && x>1);
p.addParamValue('der', 'matlab', @(x)ischar(x) && ...
    any(strcmpi(x,{'vista','matlab','kovesi'})));
p.addParamValue('hyst', [], ...
    @(x)isvector(x) && length(x)<=2 && all(x>=0) && all(x<1));
p.addParamValue('reduce', false, @(x)islogical(x) || ...
    (ischar(x) && any(strcmpi(x,{'igray','imax','isum','gmax','eor'}))));

% parse and validate all input arguments
p.parse(sig2,varargin{:}); 
p = getvarParser(p);
 
%% 
% internal variables 

C = size(I,3);

% preliminary test on the existence of the function called for edge detection
if (strcmpi(p.der,'vista') && ~exist('cannyedges','file')) || ... 
    (strcmpi(p.der,'kovesi') && (~exist('nonmaxsup','file') || ... 
    ~exist('hysthresh','file')))                                      
   error('cannyedgefuse:errormethod', ...
       ['function associated to ' p.der ' approach not found']);
end

% parameters hyst can be provided only with the function edge
if ~isempty(p.hyst) && strcmpi(p.der,'vista')
    warning('cannyedgefuse:inputwarning',['method ''matlab'' selected']); 
    p.der = 'matlab';

elseif strcmpi(p.der,'kovesi') && (isempty(p.hyst) || length(p.hyst)<2)
     error('cannyedgefuse:errorinput', ...
       ['a vector with [low,high] hysteresis values is required' ...
         ' with ' p.der ' approach']);  
end

if C==1 && ~(islogical(p.reduce) && ~p.reduce), 
    warning('cannyedgefuse:inputwarning', ...
        'ignored option ''reduce'' with scalar image');
    p.reduce = false;  
 
elseif C~=3 && ischar(p.reduce) && strcmpi(p.reduce,'igray')
    error('cannyedgefuse:inputerror', ...
         ['option ''reduce'' set to ' p.reduce ' for RGB images only']);
end

if strcmpi(p.der,'edge') && ischar(p.reduce) && strcmpi(p.reduce,'gmax')
        warning('cannyedgefuse:inputwarning', ...
            'option ''gmax'' ignored with method ''edge''');
        p.reduce = false;
end

%% 
% setting variables: check that the sigma variables are set properly, ie. 
% that |sigma1<sigma2|, otherwise invert their roles for furhter processing
if p.sig1>p.sig2                                                          
    tmp = p.sig1; p.sig1 = p.sig2; p.sig2 = tmp;
end

%% 
% main computation

edgemap = cannyedgefuse_base(I, p.sig2, p.sig1, p.serad, p.der, p.hyst, p.reduce);

%%
% display

if p.disp
    figure, imagesc(edgemap), axis image off;
    if size(edgemap,3)==1, colormap gray;  end;
    title(['fused edge map at scales ' num2str([p.sig1 p.sig2])]);
end

end % end of cannyedgefuse
