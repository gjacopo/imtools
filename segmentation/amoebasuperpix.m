%% AMOEBASUPERPIX - Amoeba-like superpixel segmentation.
%
%% Description 
% Compute superpixel regions following the approach of [GP12a]. This approach
% is similar to the Simple Linear Iterative Clustering superpixel segmentation
% technique of [ASSLFS10,LSALF10] and introducing a distance derived from the
% amoeba construction in [LDM0709]. It performs a kmeans-based clustering 
% for a selected number of seeds that outputs superpixel regions as amoeba
% neighbourhoods [LDM0709,GS08]. 
%
%% Syntax
%    Q = AMOEBASUPERPIX(I);
%    [Q, Ck, ColCk] = AMOEBASUPERPIX(I, K0, 'Property', propertyvalue, ... );
%
%% Inputs
% *|I|* : input color image of size |(X,Y,C)| (multispectral with |C=3| bands).
%
% *|K0|* : desired number of clusters, ie. of approximately equally-sized
%     superpixels.
% 
%% Property [propertyname  propertyvalues]
% *|'thres'|* : (optional) stopping criterion; it is defined as a threshold on
%     the error for relocating all the superpixel regions' centers; default:
%     |thres=eps|.
%
% *|'alpha'|* : (optional) vector of dimension |(1,n)| with |n|=1 or 2 providing
%     the influence factors used when propagating the amoeba superpixels; it
%     quantifies the relative influences of the spectral proximity w.r.t the
%     spatial one (simply computed as the local Euclidean distance) in the  
%     calculation of the perceptual distance measure; |alpha(1)| is a scalar
%     used to amplify (multiply) the spectral (graylevel or color) influence 
%     (the higher it is, the more spectral closeness is emphasized), |alpha(2)|
%     is a scalar used to introduce a factor based on the gradient magnitude 
%     (the higher it is, the less probable amoeba superpixels will cross 
%     regions of high gradient); note that the higher |alpha(1)| and/or |alpha(2)|,
%     the less compact the clusters are; default: |alpha=[1 0]|. 
%
% *|'n'|* : (optional) multiplying factor used for the search area around each
%     superpixel region's center; default: |n=2|.
%
% *|'k'|* : (optional) size of the neighbourood considered when correcting the
%     location of the superpixel regions' centers; when set to 0, the first
%     initial initialization is kept as it is; default: |k=3|.
%
% *|'iter'|* : (optional) maximum number of iterations; default: |iter=Inf|,
%     ie. the segmentation process is iterated till convergence.
%
% *|'lab'|* : (optional) boolean flag set to true for transforming, prior to
%     the processing, a 3D image (assumed to be RGB) into the Lab color
%     space; default: |lab=false|.
%
%% Outputs
% *|Q|* : diagram of amoeba superpixels regions; it takes values in the range
%     |[1,K]| where |K| is the final number of superpixel regions; it is a
%     matrix of size |(X,Y)|.
%
% *|Ck|* : coordinates of the centers of the corresponding superpixels; it is
%     a matrix of size |(K,2)|.
%
% *|ColCk|* : representative color values of the superpixels; it is a matrix
%     of size |(K,3)|.
% 
%% References
% [LDM0709]  R. Lerallut, E. Decenciere, and F. Meyer: "Image filtering using
%      morphological amoebas", Image and Vision Computing, 25(4):395-404, 2007.
%      <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.106.1413>
%
% [GS08]  J. Grazzini and P. Soille: "Adaptive morphological filters using 
%      similarities based on geodesic time", Proc. DGCI, LNCS, vol. 4992, pp.
%      519-528, 2008.
%      <http://www.springerlink.com/content/f6v62233xqkklq72/>
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
%      structures in EM images", pp. 463-471, Proc. MICCAI, 2010. 
%      <http://cvlab.epfl.ch/publications/publications/2010/LucchiSALF10.pdf>
%
%% See also
% Related:
% <SLICSUPERPIX.html |SLICSUPERPIX|>,
% <GEOSUPERPIX.html |GEOSUPERPIX|>,
% <../../propagation/html/FMM.html |FMM|>.
% Called:
% <AMOEBASUPERPIX_BASE.html |AMOEBASUPERPIX_BASE|>.

%% Function implementation
function [Q, varargout] = amoebasuperpix(I, varargin )

%%
% parsing parameters
narginchk(1, 23);
nargoutchk(1, 3);

% mandatory parameter
if ~isnumeric(I)
    error('amoebasuperpix:inputerror','a matrix is required in input'); 
end

% optional parameters
p = createParser('AMOEBASUPERPIX');   % create an instance of the inputParser class.
p.addOptional('K', 0.1, @(x)isscalar(x) && ((x>0 && x<1) || (x>1 && x==round(x))));
% additional optional parameters
p.addParamValue('alpha', [1 0], @(x) all(x>=0) && (isscalar(x) || ...
    (isnumeric(x) && length(x)<=2)));
p.addParamValue('thres', eps, @(x) isscalar(x) && x>=0);
p.addParamValue('n', 2, @(x)isscalar(x) && x>=1 && x==round(x));
p.addParamValue('k', 0, @(x)isscalar(x) && (x==0 || x>=3));
p.addParamValue('iter', Inf, @(x)isscalar(x) && x>=1);
p.addParamValue('lab', false, @islogical);
% additional parameters used by IM2POTENTIAL_BASE

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p); 

%% 
% setting variables

if size(I,3)~=3 && p.lab
    warning('amoebasuperpix:inputwarning', ...
        'option lab available for RGB and 3 vectorial image only');
    p.lab = false;
end

if p.K<1
    p.K = p.K * min(size(I,1),size(I,2)); % proportion
end

% fill p.alpha with 0's
if length(p.alpha)<2,  p.alpha = [p.alpha 0];  end

%%
% main computation
    
[Q, Ck, ColCk] = amoebasuperpix_base(I, p.K, ...
    p.lab, p.alpha, p.thres, p.n, p.k, p.iter);

if any(~Q(:))
    i = find(~Q);
    warning('amoebasuperpix:outputwarning', ...
        ['not all pixels reached by the classifier - ' ...
        'increase the range of the exploration domain (''n'') ' ...
        'or, when possible, the number of iterations (''iter'')']);
    Q(i) = 1;
end

if nargout>=2
    varargout{1} = Ck;
    if nargout>=3,  varargout{2} = ColCk;  end
end

%%
% display

if p.disp
    figure;
    if isempty(ver('images'))
        subplot(1,2,1), imagesc(rescale(I));
        subplot(1,2,2), imagesc(Q), colormap jet;
    else
        M = (imdilate(Q,ones(3,3))-Q==0);
        subplot(1,2,2), imagesc(label2rgb(Q.*M)), axis image off;
        M = cat(3,M,M,M);
        subplot(1,2,1), imagesc(rescale(I.*M)+(1-M)), axis image off;
    end
    suptitle('amoeba superpixel regions');
    figure;
    subplot(1,2,1), imagesc(rescale(I)), axis image off;
    if p.lab,  ColCk = Lab2RGB(ColCk(:,1),ColCk(:,2),ColCk(:,3));  end
    ColCk = double(reshape(ColCk(Q(:),:), size(I)));
    if exist('i','var'),  ColCk([i; i+numel(Q); i+2*numel(Q)]) = 0;  end
    subplot(1,2,2), imagesc(rescale(ColCk)), axis image off;
    suptitle('amoeba superpixel approximation');
end

end % end of amoebasuperpix
