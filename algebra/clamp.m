%% CLAMP - Clamp data.
%
%% Description
% Clamp the entry(ies) of a scalar, cell or array signal.
%
%% Syntax
%     y = CLAMP(x);
%     y = CLAMP(x, a, b);
%
%% Inputs
% *|x|* : input scalar, matrix or cell to clamp.
%
% *|a, b|* : (optional) clamping values; |a| and |b| must be scalar values 
%     or matrices of the same dimension as |x| (or |x{1}| in the case |x| is 
%     a cell); default: |a=0, b=1|.
%
%% Outputs
% *|y|* : clamped output of the same class as |x|.
%
%% See also
% Related:
% <RESCALE.html |RESCALE|>.

%% Function implementation
function y = clamp(x, a, b)

if nargin<3,  b = 1;
    if nargin<2,  a = 0;  end
end

if iscell(x),
    y = x;
    for i=1:length(x)
        y{i} = clamp(x{i}, a, b);
    end
    return;
end

[m,n] = size(x);

if ~(isequal(size(x),size(a)) || numel(a)==1) || ...
        ~(isequal(size(x),size(b)) || numel(b)==1)
    error('clamp:inputerror', ...
        'clamping value must be scalar or same size as input matrix');
end

if numel(a)==1,  a = repmat(a, [m n]);  end
if numel(b)==1,  b = repmat(b, [m n]);  end

y = max([x(:), a(:)], [], 2);
y = min([y, b(:)], [], 2);

y = reshape(y, [m,n]);
end % end of clamp
