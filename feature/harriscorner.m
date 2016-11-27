%% HARRISCORNER - Corner extraction.
%
%% Description
% Extract keypoints using Harris algorithm (with an improved version).
% 
%% Syntax
%     pt = HARRISCORNER(I);
%     pt = HARRISCORNER(I, sig, rho);
%     [pt, map] = HARRISCORNER(I, sig, rho, 'Property', propertyvalue, ...);
%     pt = HARRISCORNER(gx, gy, []);
%     [pt, map] = HARRISCORNER(gx, gy, rho, 'Property', propertyvalue, ...);
%     pt = HARRISCORNER(gx2, gy2, gxy);
%     [pt, map] = HARRISCORNER(gx2, gy2, gxy, 'Property', propertyvalue, ...);
%
%% Inputs
% Harris detector can be applied directly on the image or on already estimated
% directional derivates. In the former case, inputs are:
%      
% *|I|* : an input image with size |(X,Y,C)|, where |C=3| when |I| is multichannel.
%      
% *|sig|* : pre-smoothing width (half-window size in pixels); this parameter  
%     sets the differentiation scale in the case the image is smoothed
%     through Gaussian filtering prior to the differentiation; typ, |sig|  
%     controls the size of the objects whose orientation has to be estimated;
%     default: |sig=1|, i.e. Gaussian regularisation is used for estimating  
%     the derivatives.
%
% *|rho|* : post-smoothing width (half-window size in pixels); this parameter
%     sets the integration scale for Gaussian averaging, that controls the 
%     size of the neighbourhood in which an orientation is dominant; it is 
%     used for averaging the partial directional derivatives; default: 
%     |rho=1|.
%      
% In the latter case, the detector is applied directly on the directional
% derivatives instead, so that no redundant calculation is performed, thus
% the inputs should be: 
%      
% *|gx, gy|* : pre-computed 1st order derivatives; the calculation of Harris
%     will start from these variables; |rho| can also be passed as an
%     additional variable for computing the 2nd order derivatives.
%      
% *|gx2, gy2, gxy|* : pre-computed 2nd order derivatives; the calculation of
%     Harris will start from these variables.
%
%% Property [propertyname  propertyvalues]
% *|'method'|* : string defining the interest point response; it is either
%     |'har'| if the feature is computed following [HS88] or |'nob'| if
%     the feature is computed following [Noble89]; default: |method='har'|.
%      
% *|'kappa'|* : optional 'kappa' parameter in the case the method |'harris'|
%     is selected; default: |kappa=0.06|.
%      
% *|'thres'|* : percentage of the maximum value of the Harris corner strength 
%     used to threshold the maxima of the detector; thres is a value in |[0,1]|;
%     default: |thres=0.001|.
%      
% *|'radius'|* : size of the neighbourhood used for local maxima detection;
%     default: |radius=3|.
%
%% Outputs
% *|pt|* : matrix |(n,2)| filled with the coordinates of the |n| detected corner
%     points.
%      
% *|cornmap|* : binary mask of corner points.
%
%% References
% [HS88]  C.G. Harris and M.J. Stephens: "A combined corner and edge 
%      detector", Proc. Vision Conference, pp 147-151, 1988.
%      <http://www.bmva.org/bmvc/1988/avc-88-023.pdf>
%      
% [Noble89]  A. Noble: "Descriptions of image surfaces", PhD thesis, p45, 
%      Oxford University 1989.
%      
% [SMB00]  C. Schmid, R. Mohrand and C. Bauckhage: "Evaluation of Interest
%      Point Detectors", International Journal of Computer Vision, 
%      37(2):151-172, 2000.
%      <http://perception.inrialpes.fr/Publications/2000/SMB00/Schmid-ijcv00.pdf>
% 
%% See also
% Related:
% <CORNER.html |CORNER|>,
% <EDGECORNER.html |EDGECORNER|>,
% <SUSANCORNER.html |SUSANCORNER|>,
% <FASTCPDA.html |FASTCPDA|>,
% <FASTCORNER.html |FASTCORNER|>,
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>.
% Called: 
% <HARRISCORNER_BASE.html |HARRISCORNER_BASE|>.
 
%% Function implementation
function [pt,varargout] = harriscorner(I,varargin)

%%
% parsing and checking parameters
error(nargchk(1, 17, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

if ~isnumeric(I)
    error('harriscorner:inputerror','a matrix is required in input'); 
end

p = createParser('HARRISCORNER');   
% optional parameters
p.addOptional('V1', 1, @(x) isempty(x) || isnumeric(x) || ...
    (isscalar(x) && isfloat(x) && x>0));
p.addOptional('V2', 1, @(x) isempty(x) || isnumeric(x) || ...
    (isscalar(x) && isfloat(x) && x>0));
p.addParamValue('method', 'har', @(x)ischar(x) && ...
    any(strcmpi(x,{'nob','har','for'})));
p.addParamValue('kappa', 0.06, @(x)isscalar(x) && isfloat(x) && x>0 && x<1);
p.addParamValue('radius', 3, @(x)isscalar(x) && round(x)==x && x>0);
p.addParamValue('thres', 0.001, @(x)isscalar(x) && x>=0 && x<=1);

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% checking/setting internal variables

[X,Y,C] = size(I);                                                 

if ~isscalar(p.V1) && ~isequal(size(p.V1),size(I))
    error('harriscorner:inputerror','input matrices with ame size required');
elseif ~isscalar(p.V2) && ~isequal(size(p.V2),size(I))
    error('harriscorner:inputerror','input matrices with ame size required');
end

if any(strcmpi(p.method,{'nob','for'})),    p.kappa = 0;   end;

%%
% main processing

pt = harriscorner_base(I, p.V1, p.V2, p.kappa, p.thres, p.radius);

%%
% display

if nargout == 2 || p.disp
    tmp = false(size(I));
    for c=1:C
        tmp(pt{c}(:,1) + X*(pt{c}(:,2)-1) + X*Y*(c-1)) = 1;
    end
    
    if nargout==2,  varargout{1}=tmp;  end;
    
    if p.disp
        figure, imagesc(tmp), axis image off, title('Harris corner map');
        if size(tmp,3) == 1, colormap gray;  end;
    end
    
end % end of harriscorner



