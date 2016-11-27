%% HESSMOOTH - Hessian Tensor of a multichannel image.
%
%% Description
% Use Di Zenzo approaches for combining the different channels and various 
% different techniques for estimating the tensor.
%
%% Algorithm
% The sequential steps for computing the Hessian are: 
% 
% # smoothing, 
% # deriving,
% # resampling,
% # possibly integrating.
%
%% Syntax
%     H = HESSMOOTH(I);
%     [H, gx, gy] = HESSMOOTH(I, rho, sig);
%     [H, gx, gy, gx2, gy2, gxy] = HESSMOOTH(I, sig, ...
%                                           'Property', propertyvalue, ...);
%
%% Inputs
% *|I|* : an input image with size |(X,Y,C)|, where |C>1| when |I| is multichannel.
%
% *|rho|* : post-smoothing width (half-window size in pixels) used for defining
%     the integration scale; default: |rho=0|.
%
% *|sigma|* : smoothing width (half-window size in pixels) used for defining
%     the derivation scale; default: |sigma=1|.
%
% *|der|* : string defining the method of pre-smoothing/differentiation used
%     for estimating the directional derivatives of the input image; it is 
%     either (see |GRDSMOOTH|): |'matlab'|, |'vista'|, |'fast'|, |'conv'|,
%     |'fleck'|, |'tap5'|, |'sob'|, |'opt'| or |'ana'|; default: |der='fast'|.
%
% *|int|* : string defining the method used for the post-smoothing of the 
%     Hessian; it is either (see |GRD2GST|): |'matlab'|, |'conv'| or |'fast'| 
%     for isotropic Gaussian smoothing, or |'ani'| for anisotropic Gaussian
%     (using hourglass shaped Gaussian kernels) along the edges; this latter
%     better captures edges anisotropy; default: |int='fast'|.
%
%% Property [propertyname  propertyvalues]
% *|'samp'|* : to perform (* |samp|) interpolation of the estimated gradient
%     to avoid aliasing (should set to 2); default: no interpolation is
%     performed, therefore |samp=1|.
%
% *|'gn'|* : optional flag (|true| or |false|) to compute the tensor using
%     unit norm gradient vectors; default: |gn=false|.
%
% *|'tn'|* : optional flag (|true| or |false|) to normalize the tensor prior 
%     to its smoothing; default: |tn=false|.
%
%% Outputs
% *|H|* : an array of size |(X,Y,2,2)| storing in |H(i,j,:,:)| the Hessian 
%      at pixel |(i,j)|.
%
%% See also
% Related:
% <GSTSMOOTH.html |GSTSMOOTH|>.
% Called:
% <HESSMOOTH_BASE.html |HESSMOOTH_BASE|>.

%% Function implementation
function [H,gx,gy,gxx,gyy,gxy] = hessmooth(I,varargin)

%% 
% parsing parameters

error(nargchk(1, 25, nargin, 'struct'));
error(nargoutchk(1, 3, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(I)
    error('hessmooth:inputerror','a matrix is required in input'); 
end

% optional parameters
p = createParser('HESSMOOTH');   % create an instance of the inputParser class.
% principal optional parameters
p.addOptional('rho', 0, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addOptional('sigma', 1, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addOptional('der', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','vista','fast','conv','diag', ...
    'tap5','sob','opt','ana'}))));
p.addOptional('int', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','conv','fast','ani'}))));
% additional optional parameters
p.addParamValue('hsize',[], @(x)isscalar(x) || lenght(x)==2);
p.addParamValue('samp', 1, @(x)isscalar(x) && round(x)==x && x>=1 && x<=5);
p.addParamValue('gn', false, @(x)islogical(x));
p.addParamValue('tn', false, @(x)islogical(x));
% used with their default values only, even if it is possible to change it
p.addParamValue('thez', 8, @(x)isscalar(x) && round(x)==x);
p.addParamValue('sigt', .4, @(x)isscalar(x) && isfloat(x) && x>0);

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

% no longer here, tested in HESSMOOTH_BASE
% % define the filter window size: it is passed in order to constrain the
% % smoothing and integration window size
% if length(p.hsize)==2,   j=2;
% else                     j=1;  end
% tmp = p.hsize;   p.hsize = cell(2,1);
% if ~isempty(tmp)
%     p.hsize{1} = tmp(1);  p.hsize{2} = tmp(j);
% end
% % when p.hsize=[], then it is set a a cell of empty matrices


%% 
% main calculation

[H,gx,gy,gxx,gyy,gxy] = hessmooth_base(I, p.rho, p.sigma, p.der, p.int, ...
    p.samp, p.hsize, p.gn, p.tn, p.thez, p.sigt);

end % end of hessmooth
