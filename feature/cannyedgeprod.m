%% CANNYEDGEPROD - Canny edge detection by scale multiplication.
%
%% Description
% Two-scale Canny edge detection obtained by multiplication of Canny edge
% responses at two different scales following the approach of [BZW05].
% 
%% Syntax
%     edgemap = CANNYEDGEPROD(I, sigma2);
%     [edgemap, mag, or] = CANNYEDGEPROD(I, sigma2, sigma1, ...
%                                       'Property', propertyvalue, ...);
% 
%% Inputs
% *|I|* : input image
%
% *|sigma2, sigma1|* : standard-deviations of, resp., the coarsest and the 
%      finest Gaussian filters used for edge detection; default: |sigma2| 
%      needs always to be passed as input and |sigma1=1|.
%
%% Property [propertyname  propertyvalues]
% *|'der'|* : string setting the name of the function used for edge detection;
%      it can be: |'vista'| (or |'matlab'|) or |'kovesi'|; see function 
%      |CANNYEDGEMAP|; default: |der='matlab'|.
%
% *|'c'|* : adjustment parameter used in [BZW05] for defining the threshold
%      on the responses, see Sec.3.6; default: |c=1|.
%
% *|'reduce'|* : logical value or string defining the way the different channels
%      of a multispectral image into the edge map (see |CANNYEDGE_BASE|); it 
%      can be |true| or |false|, or one among |'igray'|, |'imax'|, |'gmax'|
%      or |'eor'|; see |CANNYEDGE_BASE| and |CANNYEDGE|.
%
%% Outputs
% *|edgemap|* : a binary edge map.
%
% *|mag, or|* : optional output images storing the magnitude and the orientation
%        of the gradient resp.
%
%% References
% [Canny86]  J.F. Canny: "A computational approach to edge detection",
%      IEEE Trans. on Pattern Analysis and Machine Intelligence, 8(6):679-698,
%      1986.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4767851&tag=1>
%
% [BZW05] P. Bao, L. Zhang and X. Wu: "Canny edge detection enhancement  
%      by scale multiplication", IEEE Trans. on Image Processing, 27(9): 
%      1485-1490, 2005. 
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1471712>
% 
%% See also
% Related:
% <CANNYEDGE.html |CANNYEDGE|>,
% <CANNYEDGEMAP.html |CANNYEDGEMAP|>,
% <EDGECORNER.html |EDGECORNER|>,
% <ROTHWELLEDGE.html |ROTHWELLEDGE|>,
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>,
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <ANISOEDGE.html |ANISOEDGE|>,
% <ELDERZUCKEREDGE.html |ELDERZUCKEREDGE|>,
% <SDGDEDGE.html |SDGDEDGE|>,
% <KOETHEDGE.html |KOETHEDGE|>,
% <PETROUEDGE.html |PETROUEDGE|>.
% Called:
% <CANNYEDGEPROD_BASE.html |CANNYEDGEPROD_BASE|>.

%% Function implementation
function [edgemap, varargout] = cannyedgeprod(I, sig2, varargin)

%% 
% parsing parameters 

error(nargchk(2, 15, nargin, 'struct'));
error(nargoutchk(1, 3, nargout, 'struct'));

if ~isnumeric(I)
    error('cannyedgeprod:inputparameter','a matrix is required in input'); 
elseif ~isscalar(sig2) || sig2<0.05
    error('cannyedgeprod:inputparameter','invalid input parameter sigma'); 
end

p = createParser('CANNYEDGEPROD');   % create an instance of the inputParser class.
p.addOptional('sig1', 1., @(x)isscalar(x) && isfloat(x) && x>0.05);
% function handling the function used for estimating Canny edges
p.addParamValue('der', 'matlab', @(x)ischar(x) && ...
    any(strcmpi(x,{'vista','matlab','kovesi'})));
p.addParamValue('c', 1, @(x)isscalar(x) && x>=0);
p.addParamValue('reduce', false, @(x)islogical(x) || ...
    (ischar(x) && any(strcmpi(x,{'igray','imax','isum','gmax','eor'}))));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% check that the sigma variables are set properly, ie. that sig1<sig2,
% otherwise invert their roles for furhter processing
if p.sig1>sig2                                                          
    tmp = p.sig1; p.sig1 = sig2; sig2 = tmp;
elseif p.sig1==sig2 
    p.sig1 =  sig2 - 0.05; 
end;

%% 
% main calculation

[edgemap, mag, or] = cannyedgeprod_base(I, sig2, p.sig1, p.der, p.c, p.reduce);

if nargout>=2, varargout{1} = mag;  end;
if nargout==3, varargout{2} = or;  end;

%%
% display

if p.disp
    figure, imagesc(edgemap), axis image off;
    if size(edgemap,3)==1, colormap gray;  end;
    title(['edge map product of scales \{' num2str(p.sig1) ':' num2str(sig2) '\}']);
end
    
end % end of cannyedgeprod
