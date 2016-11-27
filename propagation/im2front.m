%% IM2FRONT - Potential-based front propagation over an image.
%
%% Description
% Perform a (multiple) front propagation from a set of seed points using a
% metric derived from the input image over its domain.
%
%% Syntax
%   D = IM2FRONT(I, start_pts);
%   D = IM2FRONT(I, start_pts, method);
%   [D, Q] = IM2FRONT(I, start_pts, method, 'Property', propertyvalue, ...);
%
%% Inputs
% *|I|* : input image of size |(X,Y,C)|, possibly multichannel with |C>1|.
%
% *|start_pts|* : array of size |(2,k)|, where |k| is the number of starting
%     points, ie. |start_pts(:,i)| are the coordinates of the |i|-th starting
%     point.
%
% *|method|* : any string accepted by |IM2POTENTIAL|, representing the
%     method used to derive the potential function from the input image;
%     it can be: |'uni'|, |'iso'|, |'isotropic'|, |'ani'|, |'anisotropic'|,
%     |'pix'|, |'pixinv'|, |'grd'|, |'grdorth'|, |'grdn'|, |'grdninv'|,
%     |'gst'|, |'gstnorm'|, |'gstcoh'|, |'gstiso'|; see |IM2POTENTIAL| for 
%     further description.
%
%% Property [propertyname  propertyvalues]
% *|'alpha'|* : (optional) exponent factor used for amplyfying/moderating the
%     strenght of the cost function.
%
% *|'rho'|* : (optional) integration scale for spatial averaging, used for 
%     the estimation of the GST when required (see function |GSTSMOOTH|); 
%     default: |rho=3|.
%
% *|'sig'|* : differentiation scale used for the estimation of the 1st order
%     image derivatives when required (see functions |GRDSMOOTH| and 
%     |GSTSMOOTH|); default: |sig=1|.
%
% *|'der'|* : string defining the method of pre-smoothing/differentiation used
%     for estimating the directional derivatives of the input image; it is 
%     either (see |GRDSMOOTH|): |'matlab'|, |'vista'|, |'fast'|, |'conv'|, 
%     |'fleck'|, |'tap5'|, |'tap7'|, |'sob'|, |'opt'| or |'ana'|; default:
%    |der='fast'|.
%
% *|'int'|* : string defining the method used for the post-smoothing of the GST;
%     it is either (see |GRD2GST|): |'matlab'|, |'conv'|, |'fast'| or |'ani'|; 
%     default: |int='fast'|.
%
% *|'samp'|* : sampling factor used for GST estimation when required (see 
%     |GRD2GST|); default: |samp=1|.
%
% *|'eign'|* : optional string defining the method used to combine multichannel
%     in GST norm estimation when required; it is either (see |GSTFEATURE|):  
%     |'abs'|, |'l1'| (or |'zen'|), |'sum'| (or |'sap'|), |'dif'| (or |'koe'|)
%     or |'ndi'; default: |eign='l1'|. 
%
%% Outputs
% *|D|* : distance map representing the propagated front from |start_pts| using 
%     a metric derived from the image.
%
%% See also
% Related:
% <FMM_BASE.html |FMM_BASE|>,
% <IM2POTENTIAL.html |IM2POTENTIAL|>,
% <POTENTIAL2FRONT.html |POTENTIAL2FRONT|>,
% <../../derive/html/GSTSMOOTH_BASE.html |GSTSMOOTH_BASE|>, 
% Called: 
% <IM2FRONT_BASE.html |IM2FRONT_BASE|>. 

%% Function implementation
function [D, Q] = im2front(I, start_pts, varargin)

%%
% parsing and checking parameters

error(nargchk(1, 16, nargin, 'struct'));
error(nargoutchk(1, 4, nargout, 'struct'));

if ~isnumeric(I) % mandatory parameter
    error('im2front:inputerror','a matrix is required in input'); 
end

p = createParser('IM2FRONT');   % create an instance of the inputParser class.
% optional parameters
p.addOptional('method', 'gstnorm', @(x)ischar(x) && ...
    any(strcmpi(x,{'iso','isotropic','ani','anisotropic', 'isoinv',...
    'uni','pix','pixinv','grd','grdorth','grdn','grdninv','gstninv', ...
    'gst','gstnorm','gstnorm1','gstnorm2','gstnorm3','gstcoh','gstiso'})));
% additional optional parameters
p.addParamValue('alpha', [1 2], @(x)isscalar(x) || (isvector(x) && length(x)==2));
p.addParamValue('rho', 3, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addParamValue('sig', 1, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addParamValue('der', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','vista','fast','conv','fleck', ...
    'tap5','tap7','sob','opt','ana'}))));
p.addParamValue('int', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','conv','fast','ani'}))));
p.addParamValue('samp', 1, @(x)isscalar(x) && round(x)==x && x>=1 && x<=5);
p.addParamValue('eign','l1',@(x)ischar(x) && ...
    any(strcmpi(x,{'abs','zen','l1','sap','sum','ndi','dif','koe'})));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% checking/setting variables

if strcmpi(p.method,'anisotropic'),    p.method = 'ani';
elseif strcmpi(p.method,'isotropic'),  p.method = 'iso';
end

%% 
% main computation

[D, Q] = im2front_base(I, start_pts, p.method, ...
    p.alpha, p.rho, p.sig, p.der, p.int, p.samp, p.eign );

%%
% display

if p.disp
    figure;
    if isempty(ver('images'))
        subplot(1,3,1), imagesc(rescale(I));
        subplot(1,3,2), imagesc(Q), colormap jet;
    else
        M = (imdilate(Q,ones(3,3))-Q==0);
        subplot(1,3,2), imagesc(label2rgb(Q.*M)), axis image off;
        M = cat(3,M,M,M);
        subplot(1,3,1), imagesc(rescale(I.*M)+(1-M)), axis image off;
    end
    subplot(1,3,3), imagesc(rescale(D)), colormap gray, axis image off;
    hold on, plot(start_pts(1,:),start_pts(2,:),'*r'), hold off
    suptitle('front propagation');
end

end % end of im2front
