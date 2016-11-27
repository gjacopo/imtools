%% AMOEBASUPERPIX_BASE - Base function for AMOEBASUPERPIX.
%
%% Syntax
%    [Q, Ck, ColCk] = AMOEBASUPERPIX_BASE(I, K, isLab, alpha, T, n, k, maxiter);
%
%% See also
% Related:
% <AMOEBASUPERPIX.html |AMOEBASUPERPIX|>,
% <SLICSUPERPIX_BASE.html |SLICSUPERPIX_BASE|>,
% <GEOSUPERPIX_BASE.html |GEOSUPERPIX_BASE|>.
% Called:
% <matlab:web(whichpath('RGB2LAB')) |RGB2LAB|>,
% <../../graph/html/DIJKADVANCED_BASE.html |DIJKADVANCED_BASE|>,
% DIJKSTRAPROPAGATION_MEX,
% <../../propagation/html/IM2POTENTIAL.html |IM2POTENTIAL|>,
% <../../propagation/html/POTENTIAL2FRONT.html |POTENTIAL2FRONT|>,
% <../../derive/html/GRDMASK_BASE.html |GRDMASK_BASE|>,
% <../../derive/html/GSTFEATURE_BASE.html |GSTFEATURE_BASE|>,
% <../../kernel/html/NEIPOSKERNEL.html |NEIPOSKERNEL|>.

%% Function implementation
function [Q, Ck, ColCk, D] = amoebasuperpix_base(I, K, isLab, alpha, T, n, k, maxiter)

%% 
% checking/setting parameters

if nargin<8,  maxiter = Inf;
    if nargin<7,  k = 0;
        if nargin<6,  n = 2;
            if nargin<5,  T = eps;
                if nargin<4,  alpha = 1;
                end
            end
        end
    end
end

if ~exist('dijkstrapropagation_mex','file')
    warning('amoebasuperpix_base:methodwarning', ...
        ['!!!very slow matlab implementation of Disjkstra-like algorithm - use ' ...
        'C implemented version instead (see function dijkstrapropagation_mex)!!!'])
    handle_dijkstra = @(G,start) dijkadvanced(G>0, G, start);;    %#ok
    % transform the input weight matrix and create an appropriate adjacency
    % matrix for DIJKADVANCED_BASE

else
    handle_dijkstra = @(G,start) ...
        dijkstrapropagation_mex(G, start-1, -1, 1.2*numel(G));         %#ok
    % use C-implemented priority queues
end
% note: using a function handle in a loop will slower the program anyway

[X,Y,C] = size(I);

if isLab
    if C~=3                                                            
        warning('amoebasuperpix_base:inputwarning','RGB image required in input');
        if C==1,  I = repmat(I, [1 1 3]);  else  I = I(:,:,1:3);  end
        C = 3;
    end
    I = RGB2Lab(I(:,:,1), I(:,:,2), I(:,:,3));
end

%%
% for an image with |N=X*Y| pixels, the approximate size of each superpixel
% is |N/K pixels|; therefore, for roughly equally sized superpixels there
% would be a superpixel center at every grid interval |S|
S = floor(sqrt(X*Y/K));
if mod(n,2)~=0,  n = n+1;  end;
win = S * n;  pad = win / 2;  win = win + 1;

% L1 and L2 distance functions
% L1 distance 
distL1 = @(v0, v1) sum(abs(v1-v0), 2);                                 %#ok
distL2sqr = @(v0, v1) sum((v1-v0).^2, 2);
distL2 = @(v0, v1) sqrt(distL2sqr(v0,v1));


%% 
% main calculation

%%
% we begin by sampling |K| regularly spaced cluster centers
Ck = initcenters(S);
%%
% since the spatial extent of any superpixel is approximately |S^2| (the area
% of a superpixel), we can safely assume that pixels that are associated
% with this cluster center lie within a |(nS,nS)| window around the superpixel
% center on the X-Y plane
nk = size(Ck,1);

%   %----------------------------------------------------------------------
    function Ck = initcenters(S)
        % [X,Y] = size(Lab(:,:,1));
        xk = round(X / S);
        if xk==0,  xk = floor(X/2);
        else       xk = floor(S * (0:xk-1) + S/2);
        end
        yk = round(Y / S);
        if yk==0,  yk = floor(X/2);
        else       yk = floor(S * (0:yk-1) + S/2);
        end
        % ColCk = reshape(I(xk,yk,:), [length(xk)*length(yk) C]);
        [xk yk] = meshgrid(xk, yk);
        Ck = [xk(:), yk(:)];
        % Ck = sub2ind([X,Y], xk(:), yk(:));
    end
%   %----------------------------------------------------------------------

%%
% we force the centers to be located on the local extrema of the images, by
% possibly moving them to the closest local extrema.
Ck = movecenters2extremum(I, Ck);

%   %----------------------------------------------------------------------
    function Ck = movecenters2extremum(I, Ck)
        [~,lextr] = maptransition(I);
        start_pts = ind2rc([X,Y], find(lextr==C));
        [~,q] = potential2front(ones(X,Y), start_pts');
        q=q+1; % in the output of potential2front, q indexing start at 0
        Ck = start_pts(q(sub2ind([X,Y], Ck(:,1), Ck(:,2))),:);
    end
%   %----------------------------------------------------------------------

%%
% and we move the centers to seed locations corresponding to the lowest
% gradient position in a |(k,k)| neighborhood.
if k>0
    Ck = movecenters2gradmin(I, Ck, k);
end

%   %----------------------------------------------------------------------
    function Ck = movecenters2gradmin(I, Ck, k)
        [gx, gy] = grdmask_base(I,'matlab','ij');
        G = sum(gx.^2+gy.^2,3);
        % I0 = I([2:X X], (1:Y));  I1 = I([1 1:(X-1)], (1:Y));
        % G = distL2sqr(I0(:),I1(:));
        % I0 = I((1:X), [2:Y Y]);  I1 = I((1:X), [1 1:(Y-1)]);
        % G = G + distL2sqr(I0(:),I1(:));
        % G = reshape(G,[X,Y]);
        
        [ind, ix, iy] = neiposkernel(k,X);
        ind = ind(:);  ix = ix(:);  iy = iy(:);
        
        cind = sub2ind([X,Y], Ck(:,1), Ck(:,2));
        ind = repmat(cind, [1 numel(ind)]) + repmat(ind', [length(cind) 1]);
        ind(ind>X*Y) = X*Y;
        ind(ind<=0) = 1;
        [~,pos] = min(G(ind), [], 2);
        
        Ck = [Ck(:,1) + ix(pos), Ck(:,2) + iy(pos)];
        % check that the Ck are not outside the image domain's limits
        Ck(:,1) = min([max([Ck(:,1), ones(nk,1)],[],2), X*ones(nk,1)], [], 2);
        Ck(:,2) = min([max([Ck(:,2), ones(nk,1)],[],2), Y*ones(nk,1)], [], 2);
    end
%   %----------------------------------------------------------------------
 
% Following initialize the associated spectral values
cind = sub2ind([X,Y], Ck(:,1), Ck(:,2));
cind = repmat(cind, [1 C]) + (X*Y) * meshgrid(0:C-1, 1:length(cind));
ColCk = reshape(I(cind), [nk C]);


%%
% pad the original image
A = padarray(I, [pad pad 0], 'both', 'replicate');
[M,N] = size(A(:,:,1));

%%
% possibly define an additional cost function when propagating the amoeba:
% a gradient magnitude map
if alpha(2)
    rho = 1;  sigma=0.7;  samp = 1;
    int = 'fast';  der = 'fast';  eign = 'zen';
    P = gstsmooth_base(A, rho, sigma, der, int, samp, ...
        [], false, false, 8, .4);
    grd = gstfeature_base(P(:,:,1,1), P(:,:,2,2), P(:,:,1,2), ...
        'norm', eign, [], []);
    grd = rescale(grd);
    clear P;
else
    % dgxy = zeros(win^2,1);
    grd = zeros(size(A)); % we want to avoid further 'if' tests
end

%%
% reshape for easier manipulations
A = reshape(A,[M*N C]);
I = reshape(I,[X*Y C]);

%%
% initial labels and distances matrices
Q0 = zeros(M,N);
D0 = Inf(M,N);

%%
% define some utility matrices

% general indices
pixindex = reshape(1:M*N,M,N);
% relative indices of the image pixels w.r.t. the paded image 
pixin = pixindex(pad+1:pad+X,pad+1:pad+Y);

% relative indices of the pixels inside the centered analyzing window w.r.t. 
% the paded image
ind2D = neiposkernel(pad, M);  ind2D = ind2D(:);  % 2D: spatial position
%indnD = repmat(ind2D,[1 C]);            % nD: spatial/spectral position

% relative indice of the centroid pixel w.r.t. the analyzing window
cwin = [pad+1; pad+1];
cwin = sub2ind([win,win],cwin(1),cwin(2)); % = (pad+1) + (pad+1)*(win);

% list of neighbour connections within the analyzing window: we compute it
% once for all
conn = 8;  
[ic,icd] = ixneighbours(true(win,win), [], conn);
% note that ic is of size close to (but <) 8*win^2
% compute the spatial distance once for all
dxy = distL2(ind2rc(win,ic,'rc'), ind2rc(win,icd,'rc'));
% the local cost of traveling from one pixel to any of its neighbours is
% the standard Euclidean distance between them

%%
% start 'looping'
niter = 1;
err = Inf(nk,1);
% % ierr = 1:nk;

while niter<=maxiter

    % % if ~any(err>T),  break;  end
    ierr = find(err>T); 
    if isempty(ierr),  break;  end
    
    %%
    % update distance and labels
    [D0,Q0] = updatedistance(ierr, D0, Q0, Ck, ColCk);
    
    Q = Q0(pixin);
    D = D0(pixin);
    
    %%
    % update centroids
    [nCk, nColCk] = updatecentroids(ierr, Q);
    
    %%
    % update error
    
    % % ierr = 1:nk;
    err = updaterror(ierr, D, Q, Ck, nCk, ColCk, nColCk);
    
    %% 
    % update amoeba
    Ck(ierr,:) = nCk;
    ColCk(ierr,:) = nColCk;
    
    niter = niter + 1;
end

%   %----------------------------------------------------------------------
    function [D, Q, err ] = updatedistance(ierr, D, Q, Ck, ColCk)      %#ok
        % in that context, A, M, N, win, pixin, ind2D, indnD, win2D, lambda
        % pad, and dxy are all 'global' variables in the sense that they are
        % either passed to the program, or computed prior to the algorithm
        % running and they are not modified throughout the algorithm; on the
        % contrary, ierr, D, Q, Ck and ColCk are.
        % err = zeros(length(ierr),1);
        for l=1:length(ierr)  % separable estimation
            ll = ierr(l);
            % ll = l; % we recompute all of them, even those unchanged
            
            if isnan(Ck(ll,1)) && isnan(Ck(ll,2))
                break
            end
            
            % retrieve the coordinates of the cluster centroid (also center
            % of the analyzing window) w.r.t. the paded image, then the indices
            % of the analyzing window (in 2D) w.r.t. the paded image
            c2D = pixin(Ck(ll,1),Ck(ll,2));
            cind2D = c2D + ind2D;
            % % retrieve the position in color space
            % cnD = c2D + (M*N) * (0:C-1);
            % cindnD = repmat(cnD,[win^2 1]) + indnD;
            
            % extract labels and distances within the current neighbourhood
            q = Q(cind2D);  d = D(cind2D);
           
            % F = A(cindnD);
            % F(cwin + wind2D) = ColCk(ll,:);
            F = A(cind2D,:);
            F(cwin,:) = ColCk(ll,:);
            dI = distL2(F(ic,:),F(icd,:));  % distL2(F(icnD),F(icdnD));
            % the cost of traveling from ic to icd is the distance between the
            % spectral values taken by the image on those points: the closer
            % their spectral values, the higher the probability they belong
            % to the same amoeba
            
            % compute the cost induced by the gradient magnitude of the image
            % (it may be null)
            dgxy = (grd(ic) + grd(icd)) / 2;
            % the cost of traveling from ic to icd is the average of their
            % gradient values: the lower the gradient at those points, the
            % higher the probability they belong to the same amoeba; high
            % gradient pixels act like barriers in the amoeba propagation
            
            % method 1
            % start = sub2ind([win,win], pad+1, pad+1);
            W = dxy + alpha(1) * dI + alpha(2) * dgxy;
            % construct a sparse adjacency matrix G
            G = sparse(ic, icd, W, win^2, win^2);
            % Ds = handle_dijkstra(G, cwin);
            Ds = dijkstrapropagation_mex(G, cwin-1, -1, 1.2*numel(G));
            
            % % method 2
            % W = reshape(dfeat(F, repmat(F(cwin+wind2D),[win^2 1])), [win win]);
            % W = im2potential(W, 'pixinv', lambda);
            % Ds = potential2front(W, [pad+1;pad+1]);
            
            % find closest
            [d, r] = min(cat(2,d(:),Ds(:)), [], 2);
            q(r==2) = ll;
            
            % update labels
            Q(cind2D) = q; 
            
            % Qll = find(Q==ll);
            % ck = floor(sum(ind2rc(M,Qll,'rc')) / length(Qll));
            % c2D = ind2rc(M,c2D,'rc');  pos = pad+1 + ck - c2D;
            % Ds = reshape(Ds,[win,win]);
            % err(l) = Ds(pos(1), pos(2));
            
            % update distances
            D(cind2D) = d;
        end
    end
%   %----------------------------------------------------------------------
    function [Ck, ColCk] = updatecentroids(ierr, Q)
        labels = cellfun(@(k) find(Q==k), num2cell(ierr(:)), 'UniformOutput', false);
        % note that we choose the color value of the (point closest to the)
        % centroid point as a representative value
        ColCk = cellfun(@(x) sum(I(x,:),1)./length(x), labels, 'UniformOutput', false);
        % remember: for a matrix X, sum(X) and sum(X,1) are the same, hence
        % we could have written:
        %   ColCk = cellfun(@(x) sum(I(x,:))./length(x), labels, 'UniformOutput', false);
        % if we sum the Lab values over the component with same label
        %   ColCk = cellfun(@(x) sum(I([x x+X*Y x+2*X*Y]),1) / length(x), ...
        %    labels, 'UniformOutput', false);  
        labels = cellfun(@(x) ind2rc(X,x,'rc'), labels, 'UniformOutput', false);
        %disp(cellfun(@length, labels));
        Ck = cellfun(@(x) floor(sum(x,1)./size(x,1)), labels, 'UniformOutput', false);
        % note: writting here: sum(x) instead of sum(x,1) may induce errors
        % in the case classes/labels reduced to one single pixel exist
        Ck = cell2mat(Ck);
        ColCk = cell2mat(ColCk);
    end
%   %----------------------------------------------------------------------
    function err = updaterror(ierr, D, Q, Ck, nCk, ColCk, nColCk)      %#ok
        % % - error type 1:
        % err = distL1(Ck(ierr,:), nCk);
        % % - error type 2:
        ddxy =  distL2(Ck(ierr,:), nCk);
        ddrange = distL2(ColCk(ierr,:), nColCk);
        err = ddxy + alpha(1) * ddrange; % + alpha(2) * ddgxy;
        % - error type 3:
        % find the centroids' neighbours within their own superpixel
        % region (connected component with same label)
        % [i,j] = ixneighbours(true(X,Y), sub2ind([X Y],nCk(ierr,1),nCk(ierr,2)), conn);
        % ii = Q(j)~=ierr;  % Q(j)~=Q(i); the centroid may have a different label 
        % i(ii) = []; j(ii) = [];
        % % spatial distance between a centroid and ist nearest neighbours
        % ii = (ind2rc(X,i,'r')~=ind2rc(X,j,'r')) & ...
        %     (ind2rc(X,i,'c')~=ind2rc(X,j,'c'));
        % W = sqrt(2)*ii + ~ii;
        % % add the spectral distance between a centroid and any of its nearest 
        % % neighbours to the already computed distance between those neighbours
        % % and the previous centroids
        % W = W + (D(j) + lambda * distL2(nColCk(Q(i),:),I(j,:)));
        % err = accumarray(i, W, [], @(x){x}, {NaN});
        % err(cellfun(@(x)isnan(x(1)),err)) = [];
        % err = cellfun(@min,err);
    end
%   %----------------------------------------------------------------------

end % end of amoebasuperpix_base
