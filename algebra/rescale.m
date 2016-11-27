%% RESCALE - Rescale data.
%
%% Description
% Rescale a cell or array signal in |[a,b]|.
%
%% Syntax
%   y = RESCALE(x);
%   [y, d] = RESCALE(x, a, b);
%
%% Inputs
% *|x|* : input matrix or cell to rescale.
%
% *|a, b|* : (optional) rescaling values; |a| and |b| must be scalar; default:
%     |a=0, b=1|.
%
%% Outputs
% *|y|* : matrix with entries rescaled in the range |[a,b]|.
%
% *|d|* : (optional) variable returning the difference |max(x(:))-min(x(:))|
%     of the entries in the input matrix |x|.
%
%% See also
% Related:
% <CLAMP.html |CLAMP|>.

%% Function implementation
function [y, d] = rescale(x, a, b)

if nargin<3,  b = 1;
    if nargin<2,  a = 0;  end
end

if iscell(x),
    y = x;   d = cell(numel(x),1);
    for i=1:length(x)
        [y{i}, d{i}] = rescale(x{i}, a, b);
    end
    return;
end

m = min(x(:));
M = max(x(:));
d = M - m;

if d<eps,  y = x;
else         y = (b-a) * (x-m) / d + a;
end
end % end of rescale
