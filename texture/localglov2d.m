%% LOCALGLOV2D - Textural features based on the Grey Level Occurrence Vector.
% 
%% Description
% Compute local textural features using the 1D histograms of univariate 
% distribution of greylevels as defined in Grey Level Occurrence Vector (GLOV)
% approach of [GTY02,GTY03].
%
%% Syntax
%     O = LOCALGLOV2D(I);
%     O = LOCALGLOV2D(I, 'Property', propertyvalue, ...);
%
%% Inputs
% *|I|* : input image with dimension |C|, possibly multispectral (|C>1| bands).
%
%% Property [propertyname  propertyvalues]
% See function |LOCALGLCM2D| for property names and property values.
% Note that |'dcar'| and |'dpol'| are irrelevant for the GLOV approach.
%   
%% Output
% *|O|* : GLOV-based features.
%
%% References
% [GTY02]  J. Grazzini, A. Turiel, and H. Yahia: "Entropy estimation and
%      multiscale processing in meteorological satellite images", Proc. of 
%      ICPR, vol. 3, pp: 764-768, 2002.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1048103>
%
% [GTY03]  J. Grazzini, A. Turiel, and H. Yahia: "Analysis and comparison
%      of functional dependencies of multiscale textural features on 
%      monospectral infrared images", Proc. of IGARSS, vol. 3, pp: 2045-2047,
%      2003.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1294334>
%
%% See also
% <LOCALGLCM2D_BASE.html |LOCALGLCM2D_BASE|>,
% <LOCALGLSDV2D.html |LOCALGLSDV2D|>.
% Called:
% <LOCALGLOV2D_BASE.html |LOCALGLOV2D_BASE|>.

%% Function implementation
function O = localglov2d(I, varargin)

%%
% parsing parameters 

if ~isnumeric(I)
    disp('a matrix is required in input'); return;
end

% Optional parameters
p = createParser('LOCALGLOV2D');   % create an instance of the inputParser class.
% mandatory (required) variables
%p.addRequired('I', @isnumeric);
p.addParamValue('res', 64, @(x)isscalar(x) && x>=0);
p.addParamValue('win', [], @(x)isscalar(x) && x>=1);
p.addParamValue('wei', 'inv', @(x)ischar(x) && ...
    any(strcmpi(x,{'gaus','ave','inv'})));
p.addParamValue('n', 'global', @(x)ischar(x) && ...
    any(strcmpi(x,{'local','global'})));
p.addParamValue('sig', 2., @(x)isscalar(x) && x>0.1 && x<10);
p.addParamValue('mask', [], @(x)isnumeric(x) && all(x>=0));
p.addParamValue('feat', [], @(x)isnumeric(x));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% checking

if strcmp(p.wei,'gaus')
    if ~isempty(p.sig) && isempty(p.win)
        p.win = 6*p.sig + 1;
    end
elseif isempty(p.win)
    p.win = 7;
end

% ensure that the window size is odd
if ~isempty(p.win) && mod(p.win,2) == 0
    p.win = p.win + 1; 
end

%% 
% main computation

O = localglov2d_base(I, p.feat, p.res, p.win, p.wei, p.n, p.sig, p.mask );

end % end of localglov2d
