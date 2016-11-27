%% ALLCOMBS - All 'crossed' combinations of elements in various input vectors.
%
%% Description
% |A = ALLCOMBS(A1,A2,A3,...,AN)| returns all combinations of the elements
% in |A1|, |A2|, ..., and |AN|. 
%
%% Syntax
%       A = ALLCOMBS(varargin);
%
%% Inputs
% *|A1, A2, A3, ..., AN|* : inputs list of |N| matrices whose entries will
%     be combined through permutations.
%
%% Output
% *|A|* : a |(P,N)| matrix in which |P| is the product of the number of 
%     elements of the |N| inputs; empty inputs yields an empty matrix |A| 
%     of size |(0,N)|. 
%
%% Acknowledgment
% Short version of |ALLCOMB| by Jos van der Geest, no arguments' checking.
% See http://www.mathworks.com/matlabcentral/fileexchange/10064-allcomb
%
%% See also
% Related:
% <matlab:webpub(whichpath('ALLCOMB')) |ALLCOMB|>,
% <matlab:webpub(whichpath('ALLVCOMBS')) |ALLVCOMBS|>,
% <matlab:webpub(whichpath('PERMS')) |PERMS|>,
% <matlab:webpub(whichpath('NCHOOSEK')) |NCHOOSEK|>.
% Called:
% <matlab:webpub(whichpath('NDGRID')) |NDGRID|>,
% <matlab:webpub(whichpath('CAT')) |CAT|>.

%% Function implementation
function A = allcombs(varargin)  

% check for empty inputs
q = ~cellfun('isempty',varargin);

if any(~q),
    A = zeros(0,nargin) ;
else    
    ii = 1:nargin; % let the first input run fastest
    args = varargin(q) ;
    if nargin==1,  A = args{1}(:) ;
    else
        % flip using ii if last column is changing fastest: ii is used as an
        % index on the left-hand for the multiple outputs by NDGRID
        [A{ii}] = ndgrid(args{ii});  
        % concatenate into one matrix, reshape into 2D and flip columns
        A = reshape(cat(nargin+1,A{:}),[],nargin) ;
    end
end

end % end of allcombs
