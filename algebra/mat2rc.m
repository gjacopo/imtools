%% MAT2RC - Convert a matrix to a row or column vector using the |COLON| operator (:).
%
%% Syntax
%       v = MAT2RC(V, flag);
%
%% Inputs
% *|V|* : matrix of size |(X,Y,C)| to be (partially) converted to a row or 
%     column vector depending on the flag |flag| (see below); if |C>1|,
%     then the C vectorial components of |V| are converted separately.
%
% *|flag|* : (optional) string specifying the dimension the input vector
%     should be converted to, either row ('row' or 'r') or column ('col' or
%     'c'); default: |flag='c'|.
%
%% Outputs
% *|v|* : output row/column vector of size |(X*Y,C)| (case |flag='c'|) or
%     |(C,X*Y)| (case |flag='r'|).
%
%% Example
%   a=rand(4,5); mat2rc(a,'r')
%   a=rand(4,5,3); mat2rc(a,'c')
%
%% See also
% Related:
% <matlab:webpub(whichpath('COLON')) |COLON|>.

%% Function implementation
function v = mat2rc(V,flag)

if nargin<2,  flag = 'c';  end
if ~any(strcmpi(flag,{'r','row','c','col'}))
    error('mat2rc:inputerror','input flag must be ''r'' or ''c'' - see help');
end

C = size(V,3);
if C>1
    v = zeros(numel(V(:,:,1)),C);
    for ic=1:C,  v(:,ic) = mat2rc(V(:,:,ic),'c');  end
    if any(strcmpi(flag,{'r','row'})),  v = v';  end
    return
end

v = V(:);
if any(strcmpi(flag,{'r','row'})),   v = v';  end
 
end % end of mat2rc

