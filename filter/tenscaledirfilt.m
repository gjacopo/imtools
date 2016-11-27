%% TENSCALEDIRFILT - Scale adaptive directional Gaussian filtering.
%
%% Description
% Perform Gaussian adaptive directional filtering along some tensor field
% using a local adaptive scale as in [PVS06].
%
%% Syntax
%     F = TENSCALEDIRFILT(I, T);
%     F = TENSCALEDIRFILT(I, T, S);
%     F = TENSCALEDIRFILT(I, T, 'Property', propertyvalue, ...);
%
%% Inputs
% *|T|* : either a vector field (matrix of size |(n,n,2)|) or a tensor field 
%     (matrix of size |(n,n,2,2)|).
%  
%% Property [propertyname  propertyvalues]
% *|aecc|* : tuning parameter >0 setting an upper-bound on the eccentricity
%     of the filtering; default: |aecc=0.5| for a maximum eccentricity of 3 
%     when the anisotropy of the tensor field is 1.
%
% *|sig1, sig2|* : when |T| is a vector (and not tensor) field, these parameters
%     parse the (locally constant) directional scales of the anisotropic
%     Gaussian kernel: |sig1| is the (fixed) scale along the elongated 
%     orientation and is greater than or equal to |sig2|; default: |sig1=8|
%     and |sig2=3| when |T| is a vector field; otherwise, note that when |T|
%     is a tensor, the same variables are estimated (locally) from the eigen-
%     values of |T|. 
%
%% Output
% *|F|* : image filtered using the tensor field given in input.
%
%% Reference
% [PVS06]  T.Q. Pham, L.J. van Vliet, and K. Schutte: "Robust fusion of 
%      irregularly sampled data using adaptive normalized convolution",  
%      EURASIP Journal on Applied Signal Processing, 2006(1:12), ID 83268,
%      2006.
%      <http://www.hindawi.com/journals/asp/2006/083268/abs/>
% 
%% See also
% Related:
% <CONVOLUTION.html |CONVOLUTION|>,
% <TENSANIFILT.html |TENSANIFILT|>,
% <ADAPTIVEFILT.html |ADAPTIVEFILT|>,
% <GEODESICFILT.html |GEODESICFILT|>,
% <MDLFILT.html |MDLFILT|>.
% Called:
% <TENSCALEDIRFILT_BASE.html |TENSCALEDIRFILT_BASE|>.

%% Function implementation
function F = tenscaledirfilt(I,T,varargin)

%%
% parsing and checking parameters
error(nargchk(2, 13, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

if ~isnumeric(I)
    error('tenscaledirfilt:errorinput','a matrix is required in input'); 
elseif size(T,4)>2
    error('tenscaledirfilt:errorinput','T must be a tensor field or a vector field');
end

p = createParser('TENSCALEDIRFILT');   
p.addOptional('S',[], @(x)isnumeric(x) && all(x(:)>0)); 
p.addParamValue('sig1',8, @(x)isscalar(x) && x>0); 
p.addParamValue('sig2',3, @(x)isscalar(x) && x>0); 
p.addParamValue('nthe', 12, @(x)isscalar(x) && round(x)==x);
p.addParamValue('nsig', 5, @(x)isscalar(x) && round(x)==x);
p.addParamValue('nani', 4, @(x)isscalar(x) && round(x)==x);
p.addParamValue('aecc', 0.5, @(x)isscalar(x) && round(x)==x);
p.addParamValue('ani', false, @(x)islogical(x));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% main computation

F = tenscaledirfilt_base(I, T, p.S, p.sig1, p.sig2, ...
    p.nsig, p.nthe, p.nani, p.aecc, p.ani);

if p.disp
    figure, imagesc(rescale(F,0,1)), axis image off;
    if size(F,3)==1, colormap gray; end;
    title('adaptively filtered image');
end

end % end of tenscaledirfilt
