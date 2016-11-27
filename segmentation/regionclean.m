%% REGIONCLEAN - Clean up a label image.
%
%% Description
% Clean up a label image of segmentation by merging regions which do not 
% verify certain shape/area criteria to those satisfying it.
%
%% Syntax
%     S = REGIONCLEAN(L);
%     [S, Ck, ColCk]  = REGIONCLEAN(L, 'Property', propertyvalue, ...);
%
%% Inputs
% *|L : a segmentation image, where each pixel is assigned the (integer) label
%     of the region it belongs to; positive integer elements of |L| correspond
%     to different regions; labels are found in the range |[0,N]|, but are not
%     all necessarly represented.             
% 
%% Property [propertyname  propertyvalues]
% *|'area'|* : (optional) threshold |area>0| upon regions' area; when |area>1|,
%     regions whose area in pixels is below thres are merged with an adjacent
%     segment; when |0<area<1|, regions whose area is less than this treshold
%     times the mean segment area are merged with an adjacent segment;
%     default: |area=0.1| and it is estimated using |REGIONPROPS|.
% 
% *|'solid'|* : (optional) threshold |0<solid<1| upon regions' solidity (defined, 
%     for each region, as the proportion of the pixels in the convex hull of
%     this region that are also in the region and it computed as the ratio
%     Area/ConvexArea; regions whose solidity is less than this treshold are
%     merged with an adjacent segment; default: |solid=1|, ie. this is not
%     taken into account, otherwise it is estimated using |REGIONPROPS|.
% 
% *|'extent'|* : (optional) threshold |0<extent<1| upon regions' extent (defined
%     as the ratio of pixels in the region to pixels in the total bounding
%     box and computed as Area/BoundingArea); regions whose extent is less
%     than this threshold are deleted/merged; default: |extent=1|, ie. this is
%     not taken into account,  otherwise it is estimated using |REGIONPROPS|.
% 
% *|'I'|* : (optional) input image over which labels were defined; it must be
%     of same (x,y)-dimensions as the input label image |L|, it can however
%     also be multichannel (RGB or Lab); it is used in the estimation of the
%     combined spatial/spectral distance between centroids; default: |I=[]|,
%     ie. the spectral component is ignored and the distance between centroids 
%     is simply the spatial one.
% 
% *|'m'|* : (optional) scalar defining the combined spatial/spectral metric used
%     to compute distance between centroids when a small region has to be
%     merged to a larger one; the greater this value, the more spatial proximity
%     is emphasized in the calculation of distances and the more compact the 
%     clusters are; this value can be in the range |[1, 20]|; default: |m=1|. 
% 
% *|'Ck'|* : (optional) coordinates of the centroids of the corresponding regions;
%     it is a matrix of size |(N,2)|, where |N| is the maximal value of labels
%     found in |L|; in particular, non represented labels in |L| must have 
%     |NaN| entries in |Ck|; default: |Ck=[]|, centroids of the given image 
%     segmentation are estimated using |REGIONPROPS|.
% 
% *|'ColCk'|* : (optional) representative color (e.g. RGB or Lab) values of the
%     segments; it is a matrix of size |(N,3)| default: |ColCk=[]|, so that
%     the representative color is not considered when estimating the centroid
%     distance.
% 
% *|'compress'|* : (optional) flag stating if the output label image should
%     be compressed or not, ie. the range of label values is reduced; in that
%     case, |NaN| entries of the output |[Ck,ColCk]| (see below) are also
%     deleted; default: |compress=false|. 
%
%% Outputs
% *|S|* : updated segmentation image with new/refined labelled regions.
% 
% *|Ck|* : coordinates of the centroids of the regions; it is a matrix of size
%     |[K 2]|, with |K| the maximum label entry found in |S|.
% 
% *|ColCk|* : representative values of the filtered labelled regions; it is a
%     matrix of size |[K 3]|.
%
%% Example
%    I = imread('autumn.tif'); 
%    L = slicsuperpix(I, 20, 'disp', true);
%    S = REGIONCLEAN(L, 'I', I, 'area', 100, 'conn', 4, 'disp', true);
%    S = REGIONCLEAN(L, 'area', 0.1, 'conn', 8, 'm', 1); % default command
%
%% Remarks
% * The Image Processing toolbox is required. See |REGIONPROPS|.
%
% * The rule for merging regions must naturally follow decreasing predicates
%    For remind, "a (logical) predicate is said to be decreasing if and only 
%    if every subset of any set satisfying this predicate also satisfies it"
%    [SG09]. For instance, the predicate 'is the region's area less than a
%    given threshold?" (herein implemented through the use of the 'area' 
%    option) is an example of such predicate.
%
% * Regions are suppressed and/or merged at the cost of some arbitrary
%    decisions: _(i)_ in the presence of regions whose smallest distance (in 
%    the combined spatial/spectral domain) to their neighbouring regions is  
%    obtained for more than one region, _(ii)_ the regrouping process (and
%    therefore its ouput) depends on the order those regions to be merged are
%    considered.
%
% * The metric used to compute the distance between labelled regions'
%    centroids is derived from that proposed in [ASSLFS10].
% 
%% References
% [SG09]  P. Soille and J. Grazzini: "Constrained connectivity and transition 
%      regions", Proc. of ISMM, LNCS 5720, pp. 59?69, Springer-Verlag, 2009.
%      <http://www.springerlink.com/content/g6h8mk8447041532/>
%
% [ASSLFS10]  R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua and S. 
%      Susstrunk: "SLIC superpixels", EPFL Technical Report no. 149300, 2010.
%      <http://infoscience.epfl.ch/record/149300/files/SLIC_Superpixels_TR_2.pdf>
%
%% See also
% Related:
% <REGIONADJACENCY.html |REGIONADJACENCY|>,
% <matlab:webpub(whichpath('REGIONPROPS')) |REGIONPROPS|>,
% <IMLABEL.html |IMLABEL|>.
% Called:
% <REGIONCLEAN_BASE.html |REGIONCLEAN_BASE|>.

%% Function implementation
function [S, varargout] = regionclean(L, varargin)

if isempty(ver('images'))
    % note that REGIONPROPS is used in REGIONCLEAN_BASE
    error('regionclean:errortoolbox', ...
        'Image Processing toolbox required for calling REGIONPROPS');
end

%%
% parsing parameters

error(nargchk(1, 27, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(L)
    error('regionadjacency:inputerror','image of integer labels required in entry'); 
end

p = createParser('REGIONCLEAN');   % create an instance of the inputParser class.
p.addParamValue('area', 0.1, @(x)isempty(x) || (isscalar(x) && x>0));
p.addParamValue('solid', [], @(x)isempty(x) || (isscalar(x) && x>0 && x<=1));
% p.addParamValue('isoper', [], @(x)isempty(x) || (isscalar(x) && x>=0 && x<=1));
% %   'isoper' : (optional) threshold 0<isoper<1 upon regions' isoperimetric
% %     quotient (defined as the ratio of Area/CircleArea, where CircleArea is
% %     the area of the circle having the same perimeter); default: isoper=1,
% %     ie. this is not taken into account.
p.addParamValue('extent', [], @(x)isempty(x) || (isscalar(x) && x>=0 && x<=1));
% p.addParamValue('eccent', [], @(x)isempty(x) || (isscalar(x) && x>=0 && x<=1));
% %   'eccent' : (optional) threshold 0<eccent<1 upon regions' eccentricity 
% %     (specifying the eccentricity of the ellipse that has the same second
% %     moments as the region and computed as the ratio of the distance between
% %     the foci of the ellipse and its major axis length);  
p.addParamValue('conn', 8, @(x)isscalar(x) && (x==4 || x==8));
p.addParamValue('I', [], @(x)isnumeric(x));
p.addParamValue('m', 1, @(x)isscalar(x) && x>=1 && x<=20);
p.addParamValue('Ck', [], @(x)isempty(x) || (isnumeric(x) && all(x(:)>=0)));
p.addParamValue('ColCk', [], @(x)isempty(x) || isnumeric(x));
p.addParamValue('compress', false, @islogical);

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p); 

%% 
% setting variables

[X,Y] = size(L);

if ~isempty(p.Ck) && ~isequal(max(unique(L(:))),size(Ck,1))
    error('regionclean:errorinput', ...
        ['the length of the centroid matrix must be equal to the maximal ' ...
        'label value found in the input segmentation image ']);

elseif ~isempty(p.ColCk) && ~isempty(p.Ck) && ...
        ~isequal(size(ColCk,1),size(Ck,1))
    error('regionclean:errorinput', ...
        'centroid matrix and representative matrix must have same length');
end

features = {''}; thresholds = [];
if ~isempty(p.area) && p.area<X*Y,    
    features = [features, 'Area'];  
    if p.area<1, 
        props = regionprops(L,'Area');
        p.area = p.area * mean(cat(1,props.Area));
    end
    thresholds = [thresholds; p.area]; 
end 
if ~isempty(p.solid) && p.solid<1,    
    features = [features, 'Solidity'];  
    thresholds = [thresholds; p.solid];   
end
% if ~isempty(p.isoper) && p.isoper<1,  
%     features = [features, 'Perimeter', 'Area'];  
%     thresholds = [thresholds; p.isoper];
% end
% if ~isempty(p.eccent) && p.eccent<1,  
%     features = [features, 'Eccentricity'];  
%     thresholds = [thresholds; p.eccent];
% end
if ~isempty(p.extent) && p.extent<1,  
    features = [features, 'Extent'];  
    thresholds = [thresholds; p.extent];
end

features(cellfun(@isempty,features)) = []; % get rid of empty strings

if isempty(features) % and isempty(thresholds)
    error('regionclean:inputerror', ...
    'at least one discriminating feature needs to be selected');
end

% features = unique(features); % avoid repeated 'Area' entry

%% 
% main calculation

[S, p.Ck, p.ColCk] = regionclean_base(L, p.Ck, p.ColCk, p.I, p.conn, p.m, ...
    features, thresholds);

if nargout>1, varargout{1} = p.Ck; 
    if nargout>2,  varargout{2} = p.ColCk;  end
end

if p.compress,    S = compressrange(S);  end

%% 
% display

if p.disp
    figure, 
    subplot(1,2,1), imagesc(label2rgb(L.*(imdilate(L,ones(3,3))-L==0)));
    axis image off, title('input regions');
    subplot(1,2,2), imagesc(label2rgb(S.*(imdilate(S,ones(3,3))-S==0)));
    axis image off, title('cleaned regions');
end

end % end of regionclean
