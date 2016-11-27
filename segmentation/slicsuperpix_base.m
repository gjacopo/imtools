%% SLICSUPERPIX_BASE - Base function for SLICSUPERPIX.
%
%% Syntax
%    [Q, Ck, LabCk] = SLICSUPERPIX_BASE(I, K, isLab, T, win, m, k, maxiter);
%
%% See also
% Related:
% <SLICSUPERPIX.html |SLICSUPERPIX|>,
% <AMOEBASUPERPIX_BASE.html |AMOEBASUPERPIX_BASE|>,
% <GEOSUPERPIX_BASE.html |GEOSUPERPIX_BASE|>.
% Called:
% <matlab:web(whichpath('RGB2LAB')) |RGB2LAB|>,
% <../../graph/html/SCOMPONENTS.html |SCOMPONENTS|>,
% <../../kernel/html/EUCLIDKERNEL.html |EUCLIDKERNEL|>,
% <../../kernel/html/NEIPOSKERNEL.html |NEIPOSKERNEL|>.

%% Function implementation
function [Q, Ck, LabCk, D] = slicsuperpix_base(I, K, isLab, T, n, m, k, maxiter)

%%
% checking/setting parameters

if nargin<7,  maxiter = Inf;
    if nargin<6,  k = 3;
        if nargin<5,  m = 10;
            if nargin<4,  n = 2;  end
        end
    end
end

[X,Y,C] = size(I);

if C~=3
    warning('slicsuperpix_base:inputwarning','RGB image required in input');
    if C==1,   I = repmat(I, [1 1 3]);
    else       I = I(:,:,1:3); 
    end
end

%%
% for an image with N=X*Y pixels, the approximate size of each superpixel is
% N/K pixels; therefore, for roughly equally sized superpixels there would
% be a superpixel center at every grid interval S
S = floor(sqrt(X*Y/K));
if mod(n,2)~=0,  n = n+1;  end;
win = S*n;
pad = win / 2;

%%
% possibly convert the image to RGB
if isLab
    Lab = RGB2Lab(I(:,:,1), I(:,:,2), I(:,:,3));
else
    Lab = I;
end

%%
% distance functions

% L1 distance: see function UPDATERRROR below 
distL1 = @(v0, v1) sum(abs(v1-v0), 2);                                 %#ok
% L2 distance: see functions UPDATERRROR and CALCLABSPACE below
distL2 = @(v0, v1) sqrt(sum((v1-v0).^2, 2)); 

%   %----------------------------------------------------------------------
    function Ds = calclabspace(p0, p1, Lab0, Lab1, S, m) % see Eq.(1) of [ASSLFS10]       
        ddxy =  distL2(p0, p1);
        ddlab = distL2(Lab0, Lab1);
        % Ds is the sum of the lab distance and the xy plane distance normalized by
        % the grid interval S; the variable m is introduced in Ds that enables to
        % control the compactness of a superpixel.
        Ds = ddlab + m * ddxy / S;
        % note: the greater the value of m, the more spatial proximity is emphasized
        % and the more compact the cluster.
    end
%   %----------------------------------------------------------------------
% function dlab = distlab(p0, p1, Lab)
% dlab =  sqrt((Lab(p0(:,1),p0(:,2),3) -  Lab(p1(:,1),p1(:,2),3)).^2 + ...
%     (Lab(p0(:,1),p0(:,2),2) - Lab(p1(:,1),p1(:,2),2)).^2 + ...
%     (Lab(p0(:,1),p0(:,2),1) - Lab(p1(:,1),p1(:,2),1)).^2);
% end
%   %----------------------------------------------------------------------

%% 
% main calculation

%%
% we begin by sampling K regularly spaced cluster centers
[Ck, LabCk] = initcenters(Lab, S);
%%
% since the spatial extent of any superpixel is approximately S^2 (the area
% of a superpixel), we can safely assume that pixels that are associated
% with this cluster center lie within a [nS × nS] window around the superpixel
% center on the xy plane
nk = size(Ck,1);

%   %----------------------------------------------------------------------
    function [Ck, LabCk] = initcenters(Lab, S)
        % [X,Y] = size(Lab(:,:,1));
        xk = round(X / S);
        if xk==0,  xk = floor(X/2);
        else       xk = floor(S * (0:xk-1) + S/2);
        end
        yk = round(Y / S);
        if yk==0,  yk = floor(X/2);
        else       yk = floor(S * (0:yk-1) + S/2);
        end
        
        LabCk = reshape(Lab(xk,yk,:), [length(xk)*length(yk) 3]);
        
        [xk yk] = meshgrid(xk, yk);
        Ck = [xk(:), yk(:)];
        % Ck = sub2ind([X,Y], xk(:), yk(:));
    end
%   %----------------------------------------------------------------------
%     function [Ck, LabCk, S] = initcenters(I, K)
%         [Ck, S] = gridblk(I, K);
%        
%         LabCk = zeros(length(Ck),C);
%         for ic=1:C
%             x = blkproc(I(:,:,ic), S, @(x)mean(x(:)));
%             LabCk(:,ic) = x(:);
%         end
%         
%     end
% %   %--------------------------------------------------------------------

%%
% we move the centers to seed locations corresponding to the lowest gradient
% position in a [k x k] neighborhood.
if k>0
    [Ck, LabCk] = movecenters(Lab, Ck, k);
end

%   %----------------------------------------------------------------------
    function [Ck, LabCk] = movecenters(Lab, Ck, k)       
        % note that in Eq.(2), the L2 norm is used for computing the
        % distance in Lab space, not the pseudo norm defined in Eq.(1)
       
        Lab0 = Lab([2:X X], (1:Y));
        Lab1 = Lab([1 1:(X-1)], (1:Y));
        % p0 = [[2:X X]', (1:Y)'];  Lab0 = Lab(p0(:,1),p0(:,2));
        % p1 = [[1 1:(X-1)]', (1:Y)']; Lab1 = Lab(p1(:,1),p1(:,2));
        % [i,j] = ind2sub([X,Y], meshgrid(p0(:,1),p0(:,2)));
        % p0 = [i(:), j(:)];
        % [i,j] = ind2sub([X,Y], meshgrid(p1(:,1),p1(:,2)));
        % p1 = [i(:), j(:)];
        G = sum((Lab0(:)-Lab1(:)).^2, 2);
        %   distL2(Lab0(:), Lab1(:)).^2; % we avoid calling sqrt then (^2)
        %   calclabspace(p0, p1, Lab0(:), Lab1(:), S, m);
        
        Lab0 = Lab((1:X), [2:Y Y]);
        Lab1 = Lab((1:X), [1 1:(Y-1)]);
        % p0 = [(1:X)' [2:Y Y]'];      Lab0 = Lab(p0(:,1),p0(:,2));
        % p1 = [(1:X)' [1 1:(Y-1)]'];  Lab1 = Lab(p1(:,1),p1(:,2));
        % [i,j] = ind2sub([X,Y], meshgrid(p0(:,1),p0(:,2)));
        % p0 = [i(:), j(:)];
        % [i,j] = ind2sub([X,Y], meshgrid(p1(:,1),p1(:,2)));
        % p1 = [i(:), j(:)];
        G = G + sum((Lab0(:)-Lab1(:)).^2, 2);
        %       distL2(Lab0(:), Lab1(:)).^2;
        %       calclabspace(p0, p1, Lab0(:), Lab1(:), S, m);
        
        G = reshape(G,[X,Y]); % see Eq.(2) of [ASSLFS10]

        [ind, ix, iy] = neiposkernel(k,X);
        ind = ind(:);  ix = ix(:);  iy = iy(:);
        
        cind = sub2ind([X,Y], Ck(:,1), Ck(:,2));
        ind = repmat(cind, [1 numel(ind)]) + repmat(ind', [length(cind) 1]);
        ind(ind>X*Y) = X*Y;
        ind(ind<=0) = 1;
        [~,pos] = min(G(ind), [], 2);
        
        Ck = [Ck(:,1) + ix(pos), Ck(:,2) + iy(pos)];
        % check that the Ck are not outside the limit of the image
        Ck(:,1) = min([max([Ck(:,1), ones(nk,1)],[],2), X*ones(nk,1)], [], 2);
        Ck(:,2) = min([max([Ck(:,2), ones(nk,1)],[],2), Y*ones(nk,1)], [], 2);
        
        cind = sub2ind([X,Y], Ck(:,1), Ck(:,2));
        cind = [cind cind+X*Y cind+2*X*Y];
        LabCk = reshape(Lab(cind), [nk 3]);
        
    end
%   %----------------------------------------------------------------------

%%
% % define the spatial distance: it is calculated once for all: WRONG!!!
dxy = euclidkernel(win + 1);
dxy = dxy(:);
% w = numel(dxy);
w = (win+1)^2;

A = padarray(Lab, [pad pad 0], 'both', 'symmetric');
[M,N] = size(A(:,:,1));                                           

%%
% initial labels and distances matrices
Q0 = zeros(M,N);
D0 = Inf(M,N);

%%
% define some utility indices

% position indices
pixindex = reshape(1:M*N,M,N);
pixin = pixindex(pad+1:pad+X,pad+1:pad+Y);

% iindex of the centered neighbour window of analysis 
[x,y] = ndgrid(-pad:pad); ind = x + M*y;

% ind = repmat(-pad:pad,[2*pad+1 1]); ind = ind' + M*ind;
ind = ind(:);

%%
% start 'looping'
niter = 1;
err = Inf(nk,1);
ierr = 1:nk; % ierr = find(err>T);
    
while niter<=maxiter
    if ~any(err>T),  break;  end
    
    %% 
    % update labels and distance
    [D0,Q0] = updatedistance(ierr, D0, Q0, Ck, LabCk );
    
    Q = Q0(pixin);
    
    %%
    % update labels' centroids
    [nCk, nLabCk] = updatecentroids(ierr, Q);
    
    %%
    % merge regions
    closeness = 10;
    [Q, i, j] = mergeregions(ierr, Q, nCk, nLabCk, closeness);
    if ~isempty(j)
        [nCk(i,:), nLabCk(i,:)] = updatecentroids(i, Q);
        nCk(j,:) = []; nLabCk(j,:) = [];
        Ck(j,:) = [];  LabCk(j,:) = [];
        nk = size(Ck,1);
        Q = compressregions(Q, unique(Q(:)));
        Q0(pixin) = Q;
    end
    
    %% 
    % update errors
    ierr = 1:nk; % we update all of them everytime
    err = updaterror(ierr, Ck, nCk, LabCk, nLabCk);
    
    %% 
    % update labels
    Ck = nCk;
    LabCk = nLabCk;
    
    niter = niter + 1;
end

%   %----------------------------------------------------------------------
    function [D,Q] = updatedistance(ierr, D, Q, Ck, LabCk)
        for l=1:length(ierr) %nk
            % "assign the best matching pixels from a (2S × 2S) neighborhood
            % around the cluster center according to the distance measure"
            ll = ierr(l);
            % ll = l; % we recompute all of them, even those unchanged
            
            % coordinates of the neighbourhood centered around the cluster
            % center
            % if any(isnan([Ck(ll,1) Ck(ll,2)])),  continue;  end
            c = pixin(Ck(ll,1),Ck(ll,2));
            c3 =  [c, c+N*M, c+2*N*M];
            cind = c + ind;
            cind3 = repmat(c3,[w 1]) + repmat(ind,[1 3]);
            
            q = Q(cind);
            d = D(cind);
            
            dlab = distL2(A(cind3), repmat(LabCk(ll,:),[size(cind3,1) 1]));
            Ds = dlab + m * dxy; % see remark in function UPDATERROR
            
            [d, r] = min(cat(2,d,Ds(:)),[],2);
            q(r==2) = ll;
            
            % each pixel in the image is associated with the nearest cluster
            % center whose (nS × nS) search area overlaps this pixel
            Q(cind) = q;  D(cind) = d;
        end
    end
%   %----------------------------------------------------------------------
    function [Ck, LabCk] = updatecentroids(ierr, Q)
        labels = cellfun(@(k) find(Q==k), num2cell(ierr(:)), 'UniformOutput', false);
        % if we sum the Lab values over the component with same label
        % LabCk = cellfun(@(x) sum(Lab([x x+X*Y x+2*X*Y]),1) / length(x), ...
        %    labels, 'UniformOutput', false);
        labels = cellfun(@(x) [ind2rc(X,x,'r'), ind2rc(X,x,'c')], labels, ...
            'UniformOutput', false);
        Ck = cellfun(@(x) floor([sum(x(:,1)), sum(x(:,2))]/length(x(:,1))), ...
            labels, 'UniformOutput', false);
        % note that instead we choose the Lab value of the (point closest to
        % the) centroid point as a representative value
        LabCk = cellfun(@(x) squeeze(Lab(floor(x(1)),floor(x(2)),:))', ...
            Ck, 'UniformOutput', false);
        Ck = cell2mat(Ck);
        LabCk = cell2mat(LabCk);
        % if we wanted to use the Image Processing toolbox instead: no loop
        % props = regionprops(Q,'centroid', 'PixelIdxList');
        % Ck = floor(cat(1, props.Centroid));  Ck = fliplr(Ck(ierr,:));
        % LabCk = 'loop on numel(props)...';
    end
%   %----------------------------------------------------------------------
    function err = updaterror(ierr, Ck, nCk, LabCk, nLabCk)
        % algorithm 1 in [ASSLFS10] : residual error = L1 distance between
        % previous centers and recomputed centers
        % err = distL1(Ck(ierr,:), nCk);
        % instead we choose:
        err = calclabspace(Ck(ierr,:), nCk, LabCk(ierr,:), nLabCk, 1 , m);
        % note that here, by setting S=1, we replace the original expression
        % Ds = ddlab + m * ddxy / S with Ds = ddlab + m * ddxy for clarity
        % as it is easier to control the compactness of the clusters using
        % m only
    end
%   %----------------------------------------------------------------------
    function [Q, i, j] = mergeregions(ierr, Q, Ck, LabCk, closeness)
        [i,j] = find(regionadjacency_base(Q, ierr, 8)); % note: it is symmetric
        i = unique(sort([i j],2),'rows'); j = i(:,2); i = i(:,1);
        dist =  calclabspace(Ck(ierr(i),:), Ck(ierr(j),:), ...
            LabCk(ierr(i),:), LabCk(ierr(j),:), 1, m);
        l = dist < m*closeness;
        % naive merging
        if any(l)
            i = ierr(i(l)); j = ierr(j(l));
            V = sparse(i, j, ones(length(i),1), max([i(:);j(:)]),  max([i(:);j(:)]));
            V = V | V';
            [~, comp] = scomponents(V, find(sum(V,2)==1));
            i = cellfun(@(x) x(1), comp);
            j = cellfun(@(x) x(2:end), comp, 'UniformOutput', false);
            for l=1:length(i),  Q(ismember(Q,j{l})) = i(l);  end
            j = cell2mat(cellfun(@(x) x(:), j, 'UniformOutput', false));
        else
            i = []; j = [];
        end
    end
%   %----------------------------------------------------------------------
    function L = compressregions(L, labels)
        if labels(1)~=0, labels = [0; labels];  end
        V = cumsum(diff(labels)-1);
        labels = labels(2:end);
        V = labels - V;
        for l=1:length(labels),  L(L==labels(l)) = V(l);  end
    end
%   %----------------------------------------------------------------------

if nargout==4,  D = D0(pixin);  end

end % end of slicsuperpix_base





