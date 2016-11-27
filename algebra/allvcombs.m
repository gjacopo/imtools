%% ALLVCOMBS - All 'crossed' combinations of elements of a single input vector.
%
%% Description
% |A = ALLVCOMBS(A,N)| returns all combinations of N elements in vector |V|.
%
%% Syntax
%       A = ALLVCOMBS(V);
%       [A, I] = ALLVCOMBS(V, N);
%
%% Inputs
% *|V|* : array of numbers, cells or strings whose entries will be combined
%     through permutations.
%
%% Output
% *|A|* : all combinations of |N| elements of the vector |V|, with dimension
%     |(numel(V).^N,N)|. 
%
% *|I|* : optional matrix of indexes so that |A = V(I)|.
%
%% Example
%   allvcombs([0 1 1 2], 3)
%
%% See also
% Related:
% <matlab:webpub(whichpath('ALLCOMB')) |ALLCOMB|>,
% <matlab:webpub(whichpath('ALLCOMBS')) |ALLCOMBS|>,
% <matlab:webpub(whichpath('PERMS')) |PERMS|>,
% <matlab:webpub(whichpath('NCHOOSEK')) |NCHOOSEK|>.
% Called:
% <matlab:webpub(whichpath('NDGRID')) |NDGRID|>,
% <matlab:webpub(whichpath('CAT')) |CAT|>.

%% Function implementation
function [A, varargout] = allvcombs(A, N)  

if nargin<2 || isempty(N),  N = 1;  end

if isempty(A) || N == 0
    A = [] ;
    varargout{1} = [] ;

elseif fix(N) ~= N || N < 1 || numel(N) ~= 1
    error('allvcombs:errorinput', ...
        'variable ''N'' should be a >=0 integer') ;
    
elseif N==1
    A = A(:).' ;
    varargout{1} = 1:numel(A) ;
    
else
    % speed depends on the number of output arguments
    if nargout<2,  A = allvcomb(A,N) ;
    else
        % indices requested
        varargout{1} = allvcomb(1:numel(A),N) ;
        A = A(varargout{1}) ;
    end
end

%   %----------------------------------------------------------------------
    function Y = allvcomb(X,N)
        if N>1 % create a list of all possible combinations of N elements
            [Y{N:-1:1}] = ndgrid(X) ;
            Y = reshape(cat(N+1,Y{:}),[],N) ;
        else % no combinations have to be made
            Y = X(:) ;
        end
    end
%   %----------------------------------------------------------------------

end % end of allvcombs


