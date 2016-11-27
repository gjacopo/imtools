%% ANISOR - Compute the orientation information derived from the anisotropic
% continous scale model of [BBW07].
%
%% Description
%
%% Syntax
%    kappa = anisor(I);
%    kappa = anisor(I, nu, rho, sig, ...
%                            'Property', propertyvalue, ...);
%
%% Inputs
% *|I|* : an input image with size |(X,Y,C)|, where |C>1| when |I| is multichannel.
%
% *|nu|* : (optional) nonegative integer which influences the propagation of
%     the structures in the normal to the gradient direction; default: nu=4.
%
% *|rho|* : post-smoothing width; this parameter sets the integration scale
%     for spatial averaging, that controls the size of the neighbourhood
%     in which an orientation is dominant; it is  used for averaging the
%     partial directional derivatives of the tensor with a Gaussian kernel;
%     if rho<0.05, then no smoothing is performed; default: rho=1.
%
% *|sig|* : pre-smoothing width; this parameter sets the differentiation
%     scale in the case the image is smoothed prior to the differentiation 
%     through Gaussian filtering; sigma typically controls the size of the
%     objects whose orientation has to be estimated; default: sigma=1,
%     i.e. Gaussian regularisation is used for estimating the derivatives.
%
%% Property [propertyname  propertyvalues]
% *|'der'|* : string defining the method of pre-smoothing/differentiation used
%     for estimating the directional derivatives of the input image; it is 
%     either (see GRDSMOOTH): 'matlab', 'vista', 'fast', 'conv', 'diag', 
%     'tap5', 'sob', 'opt' or 'ana'; default: der='fast'.
%
% *|'int'|* : string defining the method used for the post-smoothing of the GST;
%     it is either (see GRD2GST): 'matlab', 'conv' or 'fast' for isotropic 
%     Gaussian intoothing, or 'ani' for anisotropic Gaussian (using hour-
%     glass shaped Gaussian kernels) along the edges; this latter better
%     captures edges anisotropy; default: int='fast'.
%
% *|'samp'|* : to perform (* samp) interpolation of the estimated gradient to
%     avoid aliasing (should set to 2); default: samp=1, ie no interpolation
%     is performed.
%
%% Output
% *|kappa|* : orientation information in [0,1].
%
%% Reference
% [BBW07]  M. Breuß, B. Burgeth and J. Weickert: "Anisotropic continuous-scale
%      morphology", Proc. IbPRIA, LNCS 4478, pp. 515-522, Springer, 2007.	
%      <http://www.springerlink.com/content/1hm264w86111m148/>
%
%% See also
% Related:
% <../../derive/html/GSTSMOOTH.html |GSTSMOOTH|>,
% <ANISORFILT.html |ANISORFILT|>.
% Called:
% <ANISOR_BASE.html |ANISOR_BASE|>.

%% Function implementation
function kappa = anisor(I,varargin)

%%
% parsing parameters

error(nargchk(1, 18, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(I)
    error('anisor:inputerror','a matrix is required in input'); 
end

% optional parameters
p = createParser('ANISOR');   
% principal optional parameters
p.addOptional('nu', 4, @(x)isscalar(x) && x>=0);
p.addOptional('rho', 3, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addOptional('sigma', 1, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addParamValue('der', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','vista','fast','conv','diag', ...
    'tap5','sob','opt','ana'}))));
p.addParamValue('int', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','conv','fast','ani'}))));
p.addParamValue('samp', 2, @(x)isscalar(x) && round(x)==x && x>=1 && x<=5);

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%% 
% calculation

kappa = anisor_base(I, p.nu, p.rho, p.sigma, p.der, p.int, p.samp);

%%
% display

if p.disp
    figure, imagesc(kappa), colormap gray, axis image off;
    title('orientation indice')
end

end % end of anisor
