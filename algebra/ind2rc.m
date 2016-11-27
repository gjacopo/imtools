%% IND2RC - Variant of IND2SUB behaving nicely with single output.
%
%% Description
% Output the x- and/or y- subscripts from linear index as does |IND2SUB|,
% but enable one output only.
%
%% Syntax
%    ij = IND2RC(siz, ind);
%    ij = IND2RC(siz, ind, flag);
% 
%% Inputs
% *|siz|* : dimensions of the input matrix; note that it can be reduced to
%     its x-dimension only.
% 
% *|ind|* : vector or matrix of linear indices.
% 
% *|flag|* (optional) string specifying the subscript to be computed; either
%     the row (|'r'|) or the colum (|'c'|) subscript alone will be returned;
%     flag can also be set to |'rc'| or |'cr'| to return both subscripts 
%     concatenated in the order indicated by the letter sequence; default:
%     |flag='rc'|.
%
%% Output:
% *|ij|* : matrix of dimension |(n,k)|, where |n| is the number of rows of 
%     the input |ind| variable, and |k| is equal to |length(flag)|:
% 
% * when |length(flag)=1|, |ij| is the list of row or column (depending on
%         |flag|) subscripts corresponding to the linear indices present 
%         in |ind| for a matrix of size |siz|,
% * when |length(flag)=2|, |ij| concatenates both the row and column 
%         subscripts and the order of concatenation is indicated by the
%         char sequence in |flag|.
% 
%% Remark
% This function is implemented to remedy the behaviour of |IND2SUB(siz,ind)|
% which, by default, returns |ind| when one output only is specified.
%
%% See also  
% Related: 
% <matlab:webpub(whichpath('IND2SUB')) |IND2SUB|>,
% <matlab:webpub(whichpath('SUB2IND')) |SUB2IND|>.
% Called: 
% <matlab:webpub(whichpath('FLOOR')) |FLOOR|>,
% <matlab:webpub(whichpath('MOD')) |MOD|>.

%% Function implementation
function ij = ind2rc(siz, ind, flag)

narginchk(1, 3);
nargoutchk(0, 1);

X = siz(1);

if nargin<3 || isempty(flag)
    flag = 'rc';
elseif ~(ischar(flag) && any(strcmpi(flag,{'row','col','r','c','rc','cr'})))
    error('ind2rc:errorinput', ...
        'flag must be any string among ''r'', ''c'', ''rc'' or  ''cr''')
end

ind2row = @(i) mod((i-1),X) + 1;
ind2col = @(i) floor((i-1) / X) + 1;

switch flag(1)
    case {'r','row'}
        ij = ind2row(ind);
        
    case {'c','col'}
        ij = ind2col(ind);
end

if length(flag)==2
    switch flag(2)     
        case 'r'
            ij = cat(2,ij,ind2row(ind));
            
        case 'c'
            ij = cat(2,ij,ind2col(ind));
    end
end

end % end of ind2rc