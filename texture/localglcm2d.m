%% LOCALGLCM2D - Textural features based on the Grey Level Co-occurrence Matrix. 
% 
%% Description
% Compute local textural features using the 2D histograms of joint distribution
% of greylevel pairs as defined in the Grey Level Co-occurrence Matrix (GLCM) 
% approach of [Hara79,HSD73].
%
%% Syntax
%     O = LOCALGLCM2D(I);
%     O = LOCALGLCM2D(I, 'Property', propertyvalue, ...);
%
%% Inputs
% *|I|* : input image with dimension |C|, possibly multispectral (|C>1| bands).
%
% Properties [propertyname  propertyvalues]:
% *|'dcar'|* : |(nd,2)| matrix list of nd displacement vectors, expressed in
%      cartesian coordinates X-Y, used for the estimation of the 2D joint
%      greylevel distributions; default: |nd=1| and |dcar=[1 1]|.
%
% *|'dpol'|* : similarly, |(nd,2)| matrix list of |nd| displacement vector 
%      expressed in polar coordinates; incompatible with previous |'dcar'|;
%      default: ignored, |'dcar'| is used.
%
% *|'win'|* : fixed size of the analyzing window used for local estimation; 
%
% *|'res'|* : number of levels of quantization of the image (ie. the resolution
%      of the histograms); if is set to 0, the resolution will be a constant
%      estimated w.r.t. to the window size; default: |res=64|. 
%
% *|'wei'|* : parameter defining the weighting scheme adopted for the estimation
%     of local neighbourhoods used in histogram weighting; it is either:     
%
% * |'ave'| for uniform averaging weighting function,
% * |'gaus'| for Gaussian weighting function: the neighbourhood for estimating
%          the histogram distribution is weightened by a Gaussian function
%          following the approach of [JC06],
% * |'inv'| for multiscale weighting function: the neighbourhood is weightened 
%          by the inverse euclidean distance to the central pixel following
%          the approach of [GTY02,GTY03]; 
%
% default: |wei='inv'|.
%
% *|'sig'|* : optional parameter defining either:
%
% * the standard deviation of the Gaussian window when |wei='gaus'|,
% * the scaling factor of the euclidean window when |wei='inv'|, depending
%          on the option |'wei'| selected above; 
%
% default (in both cases): |sig=2|.  
%
% *|'mask'|*: optional label image of categorisation (with values >=0).
%
% *|'im'|* : optional label(s) for which the textural features are computed
%      when the image of labels |mask| has been passed as well; |default=-1|,
%      ie. the textural features are estimated for all pixels in the input
%      image.
%
%% Output
% *|O|* : GLCM-based features.
%
%% References
% [HSD73] R.M. Haralick, K. Shanmugam, and I. Dinstein: "Textural features
%      for image classification", IEEE Trans. Systems, Man and Cybernetics,
%      3(6):610-621, 1973.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4309314>
%
% [Hara79] R.M. Haralick: "Statistical and structural approaches to texture",
%      Proceedings of IEEE, 67:786-804, 1979. 
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1455597>
%
% [ST99] L. Soh and C. Tsatsoulis: "Texture analysis of SAR sea ice imagery
%      using gray level co-occurrence matrices", IEEE Trans. on Geoscience
%      and Remote Sensing, 37(2), 1999.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=752194>
%
% [Clausi02] D A. Clausi: "An analysis of co-occurrence texture statistics 
%      as a function of grey level quantization", Canadian Journal of Remote
%      Sensing, 28(1):45-62, 2002.
%      <http://www.eng.uwaterloo.ca/~dclausi/Papers/Published%202002/Clausi%20-%20GLCP%20and%20quantization%20-%20CJRS%202002.pdf>
%
% [GTY02] J. Grazzini, A. Turiel, and H. Yahia: "Entropy estimation and
%      multiscale processing in meteorological satellite images", Proc. of 
%      ICPR, vol. 3, pp: 764-768, 2002.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1048103>
%
% [GTY03] J. Grazzini, A. Turiel, and H. Yahia: "Analysis and comparison
%      of functional dependencies of multiscale textural features on 
%      monospectral infrared images", Proc. of IGARSS, vol. 3, pp: 2045-2047,
%      2003.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1294334>
%
% [JC06] R. Jobanputra and D.A. Clausi: "Preserving boundaries for image
%      texture segmentation using grey level co-occurring probabilities",
%      Pattern Recognition, 39:234-245, 2006. 
%      <http://www.sciencedirect.com/science/article/pii/S0031320305002992>
%
%% See also
% Related:
% <LOCALGLOV2D.html |LOCALGLOV2D|>,
% <LOCALGLSDV2D.html |LOCALGLSDV2D|>.
% Called:
% <LOCALGLCM2D_BASE.html |LOCALGLCM2D_BASE|>,

%% Function implementation
function O = localglcm2d(I, varargin)

%%
% parsing parameters 
if ~isnumeric(I)
    disp('a matrix is required in input'); return;
end

% optional parameters
p = createParser('LOCALGLCM2D');   % create an instance of the inputParser class.
% mandatory (required) variables
%p.addRequired('I', @isnumeric);
p.addParamValue('res', 64, @(x)isscalar(x) && x>=0);
p.addParamValue('dcar', [1 1], @(x)isnumeric(x) && min(size(x))<=2 && ....
    all(round(x(:))==x(:)));
p.addParamValue('dpol', [], @(x)isnumeric(x) && any(size(x))==2 && ...
    all(x(:,1))>0 && -pi<all(x(:,2)) && pi>all(x(2)));
p.addParamValue('win', [], @(x)isscalar(x) && x>=1);
p.addParamValue('wei', 'ave', @(x)ischar(x) && ...
    any(strcmpi(x,{'gaus','ave','inv'})));
%p.addParamValue('n', 'global', @(x)ischar(x) && ...
%    any(strcmpi(x,{'local','global'})));
p.addParamValue('sig', '2.0', @(x)isscalar(x) && x>0.1 && x<10);
p.addParamValue('mask', [], @(x)isnumeric(x) && all(x>=0));
p.addParamValue('feat', 'all', @(x)ischar(x) && ...
    any(strcmpi(x,{'all','con','ene','ent','ent10',...
    'hom','var','dis','cor','inv'})));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% checking parameters and setting internal variables

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
% main calculation

O = localglcm2d_base(I, p.feat, p.res, p.dcar, p.dpol, p.win, p.wei, p.sig, p.mask);

end % end of localglcm2d
