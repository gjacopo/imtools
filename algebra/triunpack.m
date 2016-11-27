%% TRIUNPACK - Unpack the non-redundant representation of a triangular or symmetric matrix.
%
%% Description
% Build the square (lower or upper) triangular or symmetric matrix whose non
% redundant entries are given in a representation vector. 
%
%% Syntax
%       u = TRIUNPACK(p);
%       u = TRIUNPACK(p, x, 'l'/'u'/'s');
%
%% Inputs
% *|p|* : unpacked vector data storing the |x(x+1)/2| (see below) distinct
%     components of a (lower or upper) triangular or symmetric matrix of size
%     |(x,x)| to be output.
%
% *|x|* : (optional) size of the output (square) matrix to rebuild; default:
%     |x| is estimated as the solution (>0) of the quadratic system
%     |x(x+1)=numel(p)|.
%
% *|flag|* : (optional) string stating if the upper ('u'), the lower ('l')
%     or both ('s') parts (thus, the output matrix will be symmetric) must 
%     be filled with the unpacked data; default: |flag='s'|.
%
%% Outputs
% *|u|* : the square matrix with its lower, upper or both triangular parts
%     (depending on the variable |flag|, see above) are filled with the packed
%     representation stored in |p|.
%
%% Example
%   a = rand(4,4);
%   b = nonzeros(triu(a));
%   triunpack(b)
%   triunpack(b, [], 'u')
%   triunpack(b, [], 'l')
%
%% See also
% Called:
% <matlab:webpub(whichpath('TRIU')) |TRIU|>,
% <matlab:webpub(whichpath('TRIL')) |TRIL|>,
% <matlab:webpub(whichpath('DIAG')) |DIAG|>,
% <matlab:webpub(whichpath('OR')) |OR|>.

%% Function implementation
%--------------------------------------------------------------------------
function u = triunpack(p, x, flag)

%%
% check/set the input variables
if ~isnumeric(p) || nb_dims(p)>1
    error('triunpack:inputerror', ...
        'input unpacked data must be a numeric vector');
elseif mod(numel(p),2)
    error('triunpack:inputerror', ...
        'input unpacked data must be of odd size');    
end

if nargin<3,  flag = 's';  end
if nargin<2 || isempty(x),  x = (-1 + sqrt(1+8*numel(p))) / 2;  end

if numel(p)~=x*(x+1)/2
    error('triunpack:inputerror', ...
        'incompatible unpacked data and output matrix size');
end

%%
% create the output
u = zeros(x, x);

%%
% create the handle function
if any(strcmpi(flag,{'u','upper','s','symm'}))
    fhandle = @triu;
elseif any(strcmpi(flag,{'l','lower'}))
    fhandle = @tril;
else
    error('triunpack:inputerror', 'unreckognized flag')
end


%%
% fill the output
i = true(x,x);
u(fhandle(i)) = p;

%% 
% complete in symmetric case
if any(strcmpi(flag,{'s','symm'}))
    u(tril(i)) = p;
end

end % end of triunpack