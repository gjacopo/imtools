%% IM2AMOEBA_BASE - Base function for IM2AMOEBA.
%
%% Syntax
%     [D, Q] = IM2AMOEBA_BASE(I, start_pts, dfeat, lambda, dspace);
%
%% See also
% Related:
% <IM2AMOEBA.html |IM2AMOEBA|>, 
% <../../graph/html/DIJKSTRA_BASE.html |DIJKSTRA_BASE|>, 
% <../../graph/html/DIJK.html |DIJK|>, 
% <POTENTIAL2FRONT.html |POTENTIAL2FRONT|>,
% <FMM_BASE.html |FMM_BASE|>.
% Called: 
% <../../graph/html/DIJKADVANCED.html |DIJKADVANCED|>, 
% <../../graph/html/IXNEIGHBOURS.html |IXNEIGHBOURS|>, 
% <matlab:webpub(whichpath('SPARSE')) |SPARSE|>.

%% Function implementation
function [D, Q] = im2amoeba_base(I, start_pts, dfeat, lambda, dspace)

%%
% define internal parameters

pad = 0;
% pad = 1;
% I = padarray(I,[1 1 0], 'replicate', 'both');

[X,Y,C] = size(I);
XY = X*Y;

Q = zeros(X,Y);  Q = Q(:);
D = Inf(X,Y);  D = D(:);

%%
% main calculation over seed points

for i=1:size(start_pts,2)

    start = sub2ind([X,Y], start_pts(1,i)+pad, start_pts(2,i)+pad);
    % end_pts = [];
    % end_pts = [1:X, X*(1:Y-2)+1, X*(2:Y-1), (1:X)+X*(Y-1)];
    % note that if used with DIJKSTRAPROPAGATION, end_pts must be a row
    % vector
    % end_pts = XY;
    
    [ic,icd] = ixneighbours(I(:,:,1), [], 8);
    % ic = unique(sort([ic icd], 2), 'rows')
    % icd = ic(:,2); ic = ic(:,1);
    nrc = numel(ic);
    
    ind = repmat(XY,[nrc, 1]) * (0:C-1);
    IC = repmat(ic, [1 C]) + ind;
    ICD = repmat(icd, [1 C]) + ind;

    W = lambda * dfeat(I(IC),I(ICD));
    [x,y] = ind2sub([X,Y],ic); [IC,ICD] = ind2sub([X,Y],icd);
    W = W + dspace([x y], [IC ICD]);
    
    % construct a sparse adjacency matrix G
    G = sparse(ic, icd, W, XY, XY);
    % % transform the input weight matrix and create an appropriate adjacency
    % % matrix for DIJKADVANCED
    % A = G>0;  % adjacency matrix: connections are set to true 
    % d = dijkadvanced(A, G, start);
    niter = 1.2 * XY^2;
    d = dijkstrapropagation_mex(G, start-1, -1, niter);
   
    [D, r] = min(cat(2, d(:), D), [], 2);
    Q(r==1) = i;
    
end

D = reshape(D, [X Y]);
% D = D(1+pad:X-pad,1+pad:Y-pad);
Q = reshape(Q, [X Y]);

end % end of im2amoeba
 


