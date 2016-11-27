%% TENSANIFILT - Anisotropic Gaussian filtering of multispectral images.
%
%% Description
% Perform Gaussian adaptive filtering of an image using an anisotropic tensor 
% defined directly from the image, e.g. the Gradient Structure Tensor (GST) 
% or the Hessian.
%
%% Syntax
%     F = TENSANIFILT(I);
%     [F, S, T] = TENSANIFILT(I, method, rho, sigma, ...
%                              'Property', propertyvalue, ...);
%
%% References
% [Weick97]  J. Weickert: "Coherence-enhancing diffusion of colour images",
%      Proc. of PRIA, vol. 1, pp. 239-244, 1997.
%      <http://www.sciencedirect.com/science/article/pii/S0262885698001024>
%
% [KMS00]  R. Kimmel, R. Malladi and N. Sochen: "Images as embedded maps 
%      and minimal surfaces: movies, color, texture, and volumetric medical
%      images", International Journal of Computer Vision, 39(2):111-129,
%      2000.
%      <http://www.springerlink.com/content/u3m2629171491109/>
%
% [MN02]  M. Middendorf and H.-H. Nagel: "Empirically convergent adaptive 
%      estimation of grayvalue structure tensors", Proc. of DAGM, LNCS 2449,
%      pp. 66-74, 2002.
%      <http://www.springerlink.com/content/fbbrynkyy8yhfd4q/>
%
% [SLS01]  A.F. Sole, A. Lopez and G. Sapiro: "Crease enhancement diffusion", 
%      Computer Vision and Image Understanding, 84(2):241-248, 2001. 
%      <http://www.sciencedirect.com/science/article/pii/S1077314201909452>
%
%% See also
% Related:
% <CONVOLUTION.html |CONVOLUTION|>,
% <TENSCALEDIRFILT.html |TENSCALEDIRFILT|>,
% <ADAPTIVEFILT.html |ADAPTIVEFILT|>,
% <GEODESICFILT.html |GEODESICFILT|>,
% <MDLFILT.html |MDLFILT|>.
% Called:
% <TENSANIFILT_BASE.html |TENSANIFILT_BASE|>.

%% Function implementation
function [F, S, T] = tensanifilt(I,varargin)

%% 
% parsing parameters
error(nargchk(1, 32, nargin, 'struct'));
error(nargoutchk(1, 3, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(I)
    error('tensanifilt:inputerror','a matrix is required in input'); 
end

% optional parameters
p = createParser('TENSANIFILT');   % create an instance of the inputParser class.
p.addOptional('method', 'mid', @(x)ischar(x) && ...
    any(strcmpi(x,{'gst','gstort','wei','weickert','kim','kimmel','sol','sole',...
    'mid','middendorf','tsc','tschumperle'})));
p.addOptional('rho', 1, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addOptional('sigma', 1, @(x)isscalar(x) && isfloat(x) && x>=0);
% p.addParamValue('scales',0.5:0.1:4, @(x)nb_dims(x)==1 && length(x)>=2); 
% additional optional parameters
p.addParamValue('der', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','vista','fast','conv','fleck', ...
    'tap5','tap7','sob','opt','ana'}))));
p.addParamValue('int', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','conv','fast','ani'}))));
p.addParamValue('samp', 1, @(x)isscalar(x) && round(x)==x && x>=1 && x<=5);
% parameters for method [Weick97] and [KMS00]
p.addParamValue('a', 0.5, @(x)isscalar(x) && x>0 && x<1);
p.addParamValue('c', 1, @(x)isscalar(x) && x>0);
% parameters for method [SLS]
p.addParamValue('alpha', 1, @(x)isscalar(x) && x>=0);
p.addParamValue('beta', 1, @(x)isscalar(x) && x>=0);
p.addParamValue('eps', 0.5, @(x)isscalar(x) && x>0 && x<1);
% parameter for method [Tschumperle]
p.addParamValue('p1', 1, @(x)isscalar(x) && x>=0);
p.addParamValue('p2', 2, @(x)isscalar(x) && x>=0);

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% main computation

[F, S, T] = tensanifilt_base(I, p.method, p.rho, p.sigma, p.der, p.int, ...
    p.samp, p.a, p.c, p.alpha, p.beta, p.eps, p.p1, p.p2);

%%
% display

if p.disp
    % Display the tensor fields. The color is proportional to the size of the
    % tensor.
    U = perform_tensor_mapping(T,+1);
    U(:,:,1) = perform_histogram_equalization(U(:,:,1), 'linear');
    T1 = perform_tensor_mapping(U,-1);
    options.sub = 5;   plot_tensor_field(T1, I, options);
    figure, imagesc(rescale(F)), axis image off, title('anisotropic filtering');
    if size(F,3)==1,  colormap gray; end 
end

end % end of tensanifilt
