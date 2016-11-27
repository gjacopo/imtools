%% GSTSMOOTH - Gradient Structure Tensor (GST) of a multichannel image.
%
%% Description
% Use various different techniques for estimating the tensor and Di Zenzo
% approach [Zenzo86] for combining the different channels. 
%
%% Algorithm
% The sequential steps for computing the GST are: 
% 
% # smoothing, 
% # deriving,
% # resampling,
% # possibly integrating (spatial averaging by smoothing again).
%
%% Syntax
%     T = GSTSMOOTH(I);
%     [T, gx, gy] = GSTSMOOTH(I, rho, sig, der, int);
%     [T, gx, gy, gx2, gy2, gxy] = ....
%          GSTSMOOTH(I, rho, sig, der, int, 'Property', propertyvalue, ...);
%
%% Inputs
% *|I|* : an input image with size |(X,Y,C)|, where |C>1| when |I| is multichannel.
%
% *|rho|* : post-smoothing width; this parameter sets the integration scale
%     for spatial averaging, that controls the size of the neighbourhood
%     in which an orientation is dominant; it is  used for averaging the
%     partial directional derivatives of the tensor with a Gaussian kernel;
%     if |rho<0.05|, then no smoothing is performed; default: |rho=1|.
%
% *|sigma|* : pre-smoothing width; this parameter sets the differentiation
%     scale in the case the image is smoothed prior to the differentiation 
%     through Gaussian filtering; sigma typically controls the size of the
%     objects whose orientation has to be estimated; default: |sigma=1|,
%     i.e. Gaussian regularisation is used for estimating the derivatives.
%
% *|der|* : string defining the method of pre-smoothing/differentiation used
%     for estimating the directional derivatives of the input image; it is 
%     either (see |GRDSMOOTH|): |'matlab'|, |'vista'|, |'fast'|, |'conv'|, 
%     |'fleck'|, |'tap5'|, |'tap7'|, |'sob'|, |'opt'| or |'ana'|; default: 
%    |der='fast'|.
%
% *|int|* : string defining the method used for the post-smoothing of the GST;
%     it is either (see |GRD2GST|): |'matlab'|, |'conv'| or |'fast'| for
%     isotropic Gaussian smoothing, or |'ani'| for anisotropic Gaussian (using
%     hourglass shaped Gaussian kernels) along the edges; this latter better
%     captures edges anisotropy; default: |int='fast'|.
%
%% Property [propertyname  propertyvalues]
% *|'hsize'|* : optional filter size; default: estimated depending on |sigma|, 
%     typically |hsize=6*sigma+1|; default: |hsize=[]|, calculated later on.
%
% *|'samp'|* : to perform (* |samp|) interpolation of the estimated gradient to
%     avoid aliasing (should set to 2); default: |samp=1|, ie no interpolation
%     is performed.
%
% *|'gn'|* : optional flag (|true| or |false|) to compute the tensor using unit
%     norm gradient vectors; default: |gn=false|.
%
% *|'tn'|* : optional flag (|true| or |false|) to normalize the tensor prior to 
%     its smoothing; default: |tn=false|.
%
% *|'thez', 'sigt'|* : parameters used by the anisotropic filtering with
%     hourglass filters (see |HOURGLASSKERNEL|); default: |thez=16|, |sigt=.4|.
%
%% Outputs
% *|T|* : an array of size |(X,Y,2,2)| storing in |T(i,j,:,:)| the GST at pixel
%     |(i,j)|, ie. the symmetric  tensor that represents the convolution (|*|)
%     of a Gaussian kernel |G(rho)| with the tensor product (|.|) of the  
%     regularised spatial gradient |grad(I)| by itself [BWBM06]:
%
%           T(i,j,:,:) = G(rho) * (grad(I) . grad(I)^T) (i,j)
%
% which can be rewritten as:
%
%                          | T11(i,j)   T12(i,j) |
%           T(i,j,:,:)  =  |                     |     
%                          | T21(i,j)   T12(i,j) |
%
% with the components (|(i,j)| indexed):
%
%           T11 = G(rho) * gx2
%           T22 = G(rho) * gy2
%           T12 = T21 = G(rho) * gxy
%
% The partial derivatives |gx2|, |gy2|, and |gxy| are computed over the
%     different channels using the approach of [Zenz86]:
%
%           Gxx = sum((G(sig)*gx)^2)
%           Gyy = sum((G(sig)*gy)^2)
%           Gxy = sum((G(sig)*(gx.gy))
%
% where the sum is taken over the |C=size(I,3)| channels of the input image
%     and the filter kernel |G(sig)| ensures the regularization of the
%     gradient.
%
%% Remark
% |axis| is forced to |'xy'| in this function, ie. the orientation of the
% derivatives is the 'natural' orientation (X:horizontal, Y:vertical) unlike
% Matlab implementation (|axis='ij'|); see |GRDSMOOTH|.
%
%% References
% [Zenzo86]  S. Di Zenzo: "A note on the gradient of a multi-image",
%      Computer Vision and Graphical Image Processing, 33:116-125, 1986.
%      <http://www.sciencedirect.com/science/article/pii/0734189X86902239>
%
% [TD02]  D. Tschumperle and R. Deriche: "Diffusion PDE's on vector-valued
%      images", IEEE Signal Processing Magazine, 19(5):16-25, 2002. 
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1028349>
%
% [Scheun03]  Scheunders, P.: "A wavelet representation of multispectral
%      images", in C.H. Chen: "Frontiers of Remote Sensing Information
%      Processing", chap. 9, pp. 197-224. World Scientific, 2003.
%      http://webhost.ua.ac.be/visielab/papers/scheun/Book_Chen.pdf
%
% [Koht03a]  U. Kothe: "Edge and junction detection with an improved 
%      structure tensor", Proc. of DAGM Symposium, LNCS 2781, pp. 25-32, 
%      Springer, 2003.
%      <http://hci.iwr.uni-heidelberg.de/Staff/ukoethe/papers/structureTensor.pdf>
%
% [Koht03b]  U. Kothe: "Integrated edge and junction detection with the 
%      boundary tensor", Proc. IEEE ICCV, 2003.
%      <http://hci.iwr.uni-heidelberg.de/Staff/ukoethe/papers/polarfilters.pdf>
%
% [BWBM06]  T. Brox, J. Weickert, B. Burgeth, and P. Mrázek: "Nonlinear 
%      structure tensors", Image and Vision Computing, 24(1): 41-55, 2006. 
%      <http://www.sciencedirect.com/science/article/pii/S0262885605001642>
%
% [PC09]  G. Peyre, and L. Cohen: "Geodesic methods for shape and surface
%      processing", in "Advances in Computational Vision and Medical Image
%      Processing: Methods and Applications", vol. 13 of "Computational
%      Methods in Applied Sciences", pp. 29-56, Springer, 2009.
%
%% See also
% Related:
% <GRDSMOOTH.html |GRDSMOOTH|>,
% <HESSMOOTH.html |HESSMOOTH|>,
% <GRD2GST.html |GRD2GST|>,
% <GSTDECOMP.html |GSTDECOMP|>.
% Called:
% <GSTSMOOTH_BASE.html |GSTSMOOTH_BASE|>.

%% Function implementation
function [T,gx,gy,gx2,gy2,gxy] = gstsmooth(I,varargin)

%%
% parsing parameters

error(nargchk(1, 25, nargin, 'struct'));
error(nargoutchk(1, 6, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(I)
    error('gstsmooth:inputerror','a matrix is required in input'); 
end

% optional parameters
p = createParser('GSTSMOOTH');   
% principal optional parameters
p.addOptional('rho', 1, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addOptional('sigma', 1, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addOptional('der', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','vista','fast','conv','fleck', ...
    'tap5','tap7','sob','opt','ana'}))));
p.addOptional('int', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','conv','fast','ani'}))));
% additional optional parameters
p.addParamValue('hsize',[], @(x)isempty(x) || isscalar(x) || length(x)==2);
p.addParamValue('samp', 1, @(x)isscalar(x) && round(x)==x && x>=1 && x<=5);
p.addParamValue('gn', false, @(x)islogical(x));
p.addParamValue('tn', false, @(x)islogical(x));
% used with their default values only, even if it is possible to change it
p.addParamValue('thez', 8, @(x)isscalar(x) && round(x)==x);
p.addParamValue('sigt', .4, @(x)isscalar(x) && isfloat(x) && x>0);

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%% 
% checking variables

% reset to default
if islogical(p.int) && p.int,  p.int = 'fast';  end 
if islogical(p.der) && p.der,  p.der = 'fast';  end 

%% 
% main computation 

C = size(I,3);

[T, gx, gy, gx2, gy2, gxy] = gstsmooth_base(I, p.rho, p.sigma, ...
    p.der, p.int, p.samp, p.hsize, p.gn, p.tn, p.thez, p.sigt);

%% 
% display

if p.disp
    figure, 
    subplot(2,3,1), imagesc(rescale(gx,0,1)), axis image off, 
    title('g_x'); if C==1, colormap gray;  end
    subplot(2,3,2), imagesc(rescale(gy,0,1)), axis image off, 
    title('g_y'); if C==1, colormap gray;  end
    subplot(2,3,4), imagesc(rescale(gx2,0,1)), axis image off, 
    title('T_{11} = g_x^2'), colormap gray
    subplot(2,3,5), imagesc(rescale(gy2,0,1)), axis image off, 
    title('T_{22} = g_y^2'), colormap gray
    subplot(2,3,6), imagesc(rescale(gxy,0,1)), axis image off, 
    title('T_{12} = T_{21} = g_x \cdot g_y'), colormap gray
end

end % end of gstsmooth
