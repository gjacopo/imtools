%% GEOSUPERPIX - Geodesic superpixel segmentation.
%
%% Description
% Compute geodesic superpixelsfollowing an approach similar to the Simple
% Linear Iterative Clustering superpixel segmentation technique of
% [ASSLFS10,LSALF10]. 
%
%% Syntax
%    Q = GEOSUPERPIX(I);
%    [Q, Ck, LabCk] = GEOSUPERPIX(I, K0, method, 'Property', propertyvalue, ... );
%
%% Inputs
% *|I|* : input color image of size |(X,Y,C)| (multispectral with |C=3| bands).
%
% *|K0|* : desired number of clusters, ie. of approximately equally-sized
%     superpixels.
%
% *|method|* : string defining the method used for computing the potential
%     (and the metric) derived from the image; see function |IM2POTENTIAL|;
%     it is either based on the image itself by setting it to |'pix'| or 
%     |'pixinv'|, or based on the gradient of the image, by setting it to 
%     |'grd'|, |'grdorth'| (or |'ani'|), |'grdn'| (or |'isoinv'|) or |'grdninv'|
%     (or |'iso'|); whenever the image is a scalar (graylevel) image; when
%     the image is multispectral (|C>1|), the GST is used to derive the 
%     potential function, which will depend on the chosen method: |'iso'| (or
%     |'gstninv'|), |'gst'|, |'gstorth'|, |'ani'|, |'gstcoh'| or |'gstiso'|.
%
%% Property [propertyname  propertyvalues]
% *|'T'|* : stopping criterion; it is defined as a threshold on the errror 
%     for relocating all the superpixel regions' centers; default: |T=eps|.
%
% *|'n'|* : multiplying factor used for the search area around each superpixel
%     region's center; default: |n=2|.
%
% *|'k'|* : size of the neighbourood considered when correcting the initial
%     location of the superpixel regions' centers; when set to 0, the first
%     initialization is kept as it is; default: |k=3|.
%
% *|'iter'|* : maximum number of iterations; default: |iter=Inf|, ie. the
%     segmentation process is iterated till convergence.
%
% *|'a'|* : exponent used for amplyfying the strenght of the cost function
%     (cases |'pix'|, |'pixinv'|, |'grdn'|, |'iso'|) or the strenght of the
%     scaling function (cases |'gstninv'|, |'gstnorm'|); default: |a=1|.
%
% *|'rho'|* : integration scale for computing the GST; default: |rho=1|.
%
% *|'sigma'|* : differentiation scale for estimating the directional
%     derivatives; default: |sigma=1|.
%
% *|'der'|* : string defining the method of pre-smoothing/differentiation
%     used for estimating the directional derivatives of the input image; it 
%     is either (see |GRDSMOOTH|): |'matlab'|, |'vista'|, |'fast'|, |'conv'|, 
%     |'fleck'|, |'tap5'|, |'tap7'|, |'sob'|, |'opt'| or |'ana'|; default: 
%     |der='fast'|.
%
% *|'int'|* : string defining the method used for the post-smoothing of the 
%     GST; it is either (see GRD2GST): |'matlab'|, |'conv'| or |'fast'| for
%     isotropic Gaussian smoothing, or |'ani'| for anisotropic Gaussian (using
%     hourglass shaped Gaussian kernels) along the edges; this latter better
%     captures edges anisotropy; default: |int='fast'|.
%
% *|'samp'|* : scalar used as a sampling rate for the gradient when estimating
%     the GST; default: |samp=1|.
%
% *|'eign'|* : in the case the tensor norm estimated from the eigenvalues 
%     (|l1| and |l2|, with |l1>l2|) is to be estimated , the string |eign|
%     defines the method used for its approximation; it is either (see
%     |GSTFEATURE|): |'l1'| (or |'zen'|), |'abs'|, |'sum'| (or |'sap'|), 
%     |'dif'| (or |'koe'|) or |'ndi'|; default: |eign='l1'|. 
%
%% Outputs
% *|Q|* : diagram of geodesic superpixels regions; it takes values in the range
%     |[1,K]| where |K| is the final number of superpixel regions; it is a
%     matrix of size |(X,Y)|.
%
% *|LabCk|* : representative Lab values of the superpixels; it is a matrix 
%     of size |(K,3)|.
%
% *|Ck|* : coordinates of the centers of the corresponding superpixels; it 
%     is a matrix of size |(K,2)|.
% 
%% References
% [PC09]  G. Peyre, and L. Cohen: "Geodesic methods for shape and surface
%      processing", in "Advances in Computational Vision and Medical Image
%      Processing: Methods and Applications", vol. 13 of "Computational 
%      Methods in Applied Sciences", pp. 29-56, Springer, 2009.
%      <http://www.springerlink.com/content/k15pkj1434003774/>
%
% [GSD10]  J. Grazzini, P. Soille and S. Dillard: "Multichannel image 
%      regularisation using anisotropic geodesic filtering", Proc. ICPR,
%      pp. 2664-2667, 2010.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5596008>
%
% [ASSLFS10]  R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua and S. 
%      Susstrunk: "SLIC superpixels", EPFL Technical Report no. 149300, 2010.
%      <http://infoscience.epfl.ch/record/149300/files/SLIC_Superpixels_TR_2.pdf>
%
% [LSALF10]  A. Lucchi, K. Smith, R. Achanta, V. Lepetit and P. Fua: "A 
%      fully automated approach to segmentation of irregularly shaped cellular
%      structures in EM images", Proc. MICCAI, 2010. 
%      <http://cvlab.epfl.ch/publications/publications/2010/LucchiSALF10.pdf>
%
%% See also
% Related:
% <SLICSUPERPIX.html |SLICSUPERPIX|>,
% <AMOEBASUPERPIX.html |AMOEBASUPERPIX|>,
% <../../filter/html/GEODESICFILT.html |GEODESICFILT|>,
% <../../propagation/html/FMM.html |FMM|>.
% Called:
% <GEOSUPERPIX_BASE.html |GEOSUPERPIX_BASE|>.

%% Function implementation
function [Q, varargout] = geosuperpix(I, varargin)

%%
% parsing parameters
error(nargchk(1, 37, nargin, 'struct'));
error(nargoutchk(1, 3, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(I)
    error('geosuperpix:inputerror','a matrix is required in input'); 
end

% optional parameters
p = createParser('GEOSUPERPIX');   % create an instance of the inputParser class.
p.addOptional('K', 0.1, @(x)isscalar(x) && ((x>0 && x<1) || (x>1 && x==round(x))));
p.addOptional('method', 'ani',  @(x)ischar(x) && ...
    any(strcmpi(x,{'iso','isotropic','ani','anisotropic', ...
    'gstninv','gst','gstorth','gstn','gstn1','gstn2','gstn3','gstcoh'})));
% additional optional parameters
p.addParamValue('T', eps, @(x) isscalar(x) && x>=0);
p.addParamValue('n', 2, @(x)isscalar(x) && x>=1 && x==round(x));
p.addParamValue('k', 3, @(x)isscalar(x) && (x==0 || x>=3));
p.addParamValue('iter', Inf, @(x)isscalar(x) && x>=1);
% additional parameters used by IM2POTENTIAL_BASE
p.addParamValue('rho', 1, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addParamValue('sigma', 1, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addParamValue('der', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','vista','fast','conv','fleck', ...
    'tap5','tap7','sob','opt','ana'}))));
p.addParamValue('int', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','conv','fast','ani'}))));
p.addParamValue('samp',1, @(x)isscalar(x) && round(x)==x && x>=1);
p.addParamValue('eign','l1',@(x)ischar(x) && ...
    any(strcmpi(x,{'abs','zen','l1','sap','sum','ndi','dif','koe'})));
p.addParamValue('a', [1 1], @(x)isscalar(x) || ...
    (isvector(x) && length(x)==2));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p); 

%%
% setting variables

if p.K<1
    p.K = p.K * min(size(I,1),size(I,2)); % proportion
end
   
%% 
% main computation

[Q, Ck, LabCk] = ...
    geosuperpix_base(I, p.K, p.T, p.n, p.k, p.iter, p.method, ...
    p.rho, p.sigma, p.a, p.der, p.int, p.samp, p.eign);

if nargout>=2
    varargout{1} = Ck;
    if nargout>=3,  varargout{2} = LabCk;  end
end

%% 
% display

if p.disp
    figure;
    if isempty(ver('images')),       imagesc(Q), colormap jet;
    else  imagesc(label2rgb(Q.*(imdilate(Q,ones(3,3))-Q==0)));
    end
    axis image off, title('geodesic superpixel areas');
    rgb = Lab2RGB(LabCk(:,1),LabCk(:,2),LabCk(:,3));
    rgb = double(reshape(rgb(Q(:),:), size(I)));
    figure, imagesc(rescale(rgb)), axis image off, 
    title('geodesic superpixel representation');
end

end % end of geosuperpix