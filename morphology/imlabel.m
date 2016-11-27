%% IMLABEL - Labelling of connected components in 2-D arrays
%
%% Description
% Label connected components in 2-D arrays. It is a generalization of |BWLABEL|:
% |BWLABEL| works with 2-D binary images only, whereas |IMLABEL| works with 2-D
% arrays of any class. 
%
%% Syntax
%     L = IMLABEL(I);
%     [L, num, sz] = IMLABEL(I, n);
%
%% Inputs
% *|I|* : input (logical or numeric) image.
%
% *|N|* : (optional) integer variable specifying the graph connectivity; it
%     can have a value of either 4 (4-connected objects) or 8 (8-connected 
%     objects); default: |N=8|. 
%    
%% Outputs
% *|L|* : matrix of the same size as |I|, containing labels for the connected
%     components in |I|; the elements of |L| are integer values greater than
%     or equal to 0; the pixels labeled 0 are the background, corresponding
%     to the |NaN| components of the input array; two adjacent components
%     (pixels), of respective indexes |IDX1| and |IDX2|, are connected if
%     |I(IDX1)| and |I(IDX2)| are equal.  
%
% *|num|* : (optional) number of connected objects found in |I|.   
%
% *|sz|* : matrix of the same size as |I| that contains the sizes of the
%     connected objects: for a pixel whose index is |IDX|, we have: 
%     |sz(IDX) = nnz(L==L(IDX))|.
%
%% Examples
%    I = [3 3 3 0 0 0 0 0
%         3 3 1 0 6.1 6.1 9 0
%         1 3 1 3 6.1 6.1 0 0
%         1 3 1 3 0 0 1 0
%         1 3 3 3 3 3 1 0
%         1 3 1 0 0 3 1 0
%         1 3 1 0 0 1 1 0
%         1 1 1 1 1 0 0 0];
%    L4 = IMLABEL(I,4);
%    L8 = IMLABEL(I,8);
%    subplot(211), imagesc(L4), axis image off
%    title('pixels of same color belong to the same region (4-connection)')
%    subplot(212), imagesc(L8), axis image off    
%    title('pixels of same color belong to the same region (8-connection)')
%    % Comparison between BWLABEL and LABEL:
%    BW = logical([1 1 1 0 0 0 0 0
%                  1 1 1 0 1 1 0 0
%                  1 1 1 0 1 1 0 0
%                  1 1 1 0 0 0 1 0
%                  1 1 1 0 0 0 1 0
%                  1 1 1 0 0 0 1 0
%                  1 1 1 0 0 1 1 0
%                  1 1 1 0 0 0 0 0]);
%    L = bwlabel(BW,4);
%    % The same result can be obtained with LABEL:
%    BW2 = double(BW);
%    BW2(~BW) = NaN;
%    L2 = label(BW2,4);
%
%% Acknowledgments
% * This code is due to Damien Garcia -- <http://www.biomecardio.com>
% * The Union-Find algorithm is based on the following document:
% <http://www.cs.duke.edu/courses/cps100e/fall09/notes/UnionFind.pdf>
%
%% Remarks
% * Use |BWLABEL| if the input is binary since |BWLABEL| will be much faster.
% * Note that |NaN| values are ignored and considered as background; as |IMLABEL|
% works with arrays of any class, the 0s are NOT considered as the background. 
%
%% See also
% Related:
% <.html ||>,
% <../..//html/.html ||>,
% <matlab:web(whichpath('')) ||>,
% Related: BWLABEL, BWLABELN, LABEL2RGB.

%% Function implementation
function [L,num,sz] = imlabel(I,n)

%% 
% check input arguments
error(nargchk(1,2,nargin));
if nargin==1, n=8; end

assert(ndims(I)==2,'The input I must be a 2-D array')

%%
% initialization of the two arrays (ID & SZ) required during the
% Union-Find algorithm.
sizI = size(I);
id = reshape(1:prod(sizI),sizI);
sz = ones(sizI);

% indexes of the adjacent pixels
vec = @(x) x(:);
if n==4 % 4-connected neighborhood
    idx1 = [vec(id(:,1:end-1)); vec(id(1:end-1,:))];
    idx2 = [vec(id(:,2:end)); vec(id(2:end,:))];
elseif n==8 % 8-connected neighborhood
    idx1 = [vec(id(:,1:end-1)); vec(id(1:end-1,:))];
    idx2 = [vec(id(:,2:end)); vec(id(2:end,:))];
    idx1 = [idx1; vec(id(1:end-1,1:end-1)); vec(id(2:end,1:end-1))];
    idx2 = [idx2; vec(id(2:end,2:end)); vec(id(1:end-1,2:end))];
else
    error('The second input argument must be either 4 or 8.')
end

%%
% create the groups and merge them (Union/Find Algorithm)
for k = 1:length(idx1)
    root1 = idx1(k);
    root2 = idx2(k);
    
    while root1~=id(root1)
        id(root1) = id(id(root1));
        root1 = id(root1);
    end
    while root2~=id(root2)
        id(root2) = id(id(root2));
        root2 = id(root2);
    end
    
    if root1==root2, continue, end
    % (The two pixels belong to the same group)
    
    N1 = sz(root1); % size of the group belonging to root1
    N2 = sz(root2); % size of the group belonging to root2
    
    if I(root1)==I(root2) % then merge the two groups
        if N1 < N2
            id(root1) = root2;
            sz(root2) = N1+N2;
        else
            id(root2) = root1;
            sz(root1) = N1+N2;
        end
    end
end

while 1
    id0 = id;
    id = id(id);
    if isequal(id0,id), break, end
end
sz = sz(id);

%%
% label matrix
isNaNI = isnan(I);
id(isNaNI) = NaN;
[id,~,n] = unique(id);
I = 1:length(id);
L = reshape(I(n),sizI);
L(isNaNI) = 0;

if nargin>1, num = nnz(~isnan(id)); end

end % end of imlabel
