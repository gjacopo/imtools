%% GEOSUPERPIX_BASE - Base function for GEOSUPERPIX.
%
%% Syntax
%    [Q, LabCk, Ck] = geosuperpix_base(I, K, T, n, k, maxiter, method, ...
%                                 rho, sigma, a, der, int, samp, eign);
%
%% See also
% Related:
% <GEOSUPERPIX.html |GEOSUPERPIX|>,
% <AMOEBASUPERPIX_BASE.html |AMOEBASUPERPIX_BASE|>,
% <SLICSUPERPIX_BASE.html |SLICSUPERPIX_BASE|>.
% Called:
% <matlab:web(whichpath('RGB2LAB')) |RGB2LAB|>,
% <../../propagation/html/IM2POTENTIAL.html |IM2POTENTIAL|>,
% <../../propagation/html/POTENTIAL2FRONT.html |POTENTIAL2FRONT|>,
% <../../derive/html/GSTDECOMP.html |GSTDECOMP|>.

%% Function implementation
%--------------------------------------------------------------------------
function [Q, Ck, LabCk] = ...
    geosuperpix_base(I, K, T, n, k, maxiter, method, ...
    rho, sigma, a, der, int, samp, eign)

% compute the set of gradient minima and their influence zone
% select seeds on an uniform grid

% for iteration
% compute the distance to the seeds
% end

% we accept a variable number of inputs
if nargin<14,     eign = 'l1';
    if nargin<13,     samp = 1;
        if nargin<12,     int = 'fast';
            if nargin<11,     der = 'fast';
                if nargin<10,     a = [1 2];
                    if nargin<9,     sigma = 1;
                        if nargin<8,     rho = 1;
                            if nargin7,    method = 'ani';
                                if nargin<6,  maxiter = Inf;
                                    if nargin<5,  k = 3;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

[X,Y,C] = size(I);                                                 

if C~=3
    warning('geosuperpix_base:inputwarning','RGB image required in input');
    if C==1,   I = repmat(I, [1 1 3]);
    else       I = I(:,:,1:3); 
    end
end

%%
% distance
distL2 = @(v0, v1) sqrt(sum((v1-v0).^2, 2));

%%
% for an image with N=X*Y pixels, the approximate size of each superpixel is
% N/K pixels; therefore, for roughly equally sized superpixels there would 
% be a superpixel center at every grid interval S:
S = floor(sqrt(X*Y/K));
if mod(n,2)~=0,  n = n+1;  end;
win = S*n;
pad = win / 2;

Lab = RGB2Lab(I(:,:,1), I(:,:,2), I(:,:,3));

A = padarray(Lab, [pad pad 0], 'both', 'replicate');
[M,N] = size(A(:,:,1));                                           

%%
% we begin by sampling K regularly spaced cluster centers
[Ck, LabCk] = initcenters(I, K);
nk = length(Ck);

% % and we move the centers to seed locations corresponding to the lowest 
% % gradient position in a [k x k] neighborhood.
% if k>0
%     [Ck, LabCk] = movecenters(distL2, Lab, Ck, k, S, m);
% end


%%
% compute once for all the potential function based on the gradient and/or
% the gradient structure tensor of the multispectral input image and define 
% the propagation function: it is either a scalar field (isotropic case) 
% providing with the speed (strenght) for the propagation (the higher the 
% value at one point, the faster the propagation through this point) or a 
% tensor field (anisotropic case) providing with a speed and a direction for
% the propagation.
     
[TT, L] = im2potential(A, method, a, rho, sigma, der, int, samp, eign); %#ok

[~, L2, E1, E2] = gstdecomp(TT);

%A = padarray(Lab, [pad pad 0], 'both', 'replicate');
%[M,N] = size(A(:,:,1));                                           

Q = zeros(X,Y);
D = Inf(M,N); % D = D(:);

% indexes
pixindex = reshape(1:M*N,M,N);
pixin = pixindex(pad+1:pad+X,pad+1:pad+Y);

% Index of the centered neighbour window of analysis 
ind = repmat(-pad:pad,[2*pad+1 1]);
[x,y] = size(ind);
ind = ind' + M*ind;
ind = ind(:);

% start looping
niter = 1;
err = Inf(nk,1);

while niter<=maxiter
    
    ierr = find(err>T);
    
    if isempty(ierr),  break;  end
         
    Q = padarray(Q, [pad pad], 'both', 'replicate');
   
    for k=1:nk %length(ierr)
 
        kk = k; 

        c = pixin(Ck(kk,1),Ck(kk,2));
        c3 =  [c, c+N*M, c+2*N*M];
        cind = c + ind;
        cind2 = repmat(c3(1:2),[x*y 1]) + repmat(ind,[1 2]);
        cind3 = repmat(c3,[x*y 1]) + repmat(ind,[1 3]);
        q = Q(cind);
        d = D(cind);
        l2 = reshape(L2(cind), [x y]);  % l1 = reshape(L1(cind), [x y]);
        e2 = reshape(E2(cind2), [x y 2]);  e1 = reshape(E1(cind2), [x y 2]);
        
        % a = reshape(distL2(A(cind3), squeeze(repmat(LabCk(k,:),[x*y 1]))), [x y]);
         a = reshape(distL2(A(cind3), repmat(LabCk(k,:),[x*y 1])), [x y]);
        % no, should average value
        l1 = 1 ./ (1+a);
        l2 = l2 ./ (1+a);

        ss = gstdecomp(l1, l2, e1, e2);
        dd = potential2front(ss, [pad+1;pad+1]);
        
        [d, r] = min(cat(2,d(:),dd(:)), [], 2);
        q(r==2) = k;
        
        Q(cind) = q;
        D(cind) = d;

    end
    
    Q = Q(pad+1:pad+X,pad+1:pad+Y);

    [nCk, nLabCk] = updatecenters(ierr, Q, Lab);
    warning('geosuperpix_base:warning','!!!update error vector!!!');
    
    Ck(ierr,:) = nCk;
    LabCk(ierr,:) = nLabCk;
    
    niter = niter + 1;
end

end % end of geosuperpix_base


%% Subfunction

%--------------------------------------------------------------------------
function [Ck, LabCk] = initcenters(Lab, S)
[X,Y] = size(Lab(:,:,1));
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
% %----------------------------------------------------------------------
% function [Ck, LabCk, S] = initcenters(I, K)
% [Ck, S] = gridblk(I, K);
% 
% LabCk = zeros(length(Ck),C);
% for ic=1:C
%     x = blkproc(I(:,:,ic), S, @(x)mean(x(:)));
%     LabCk(:,ic) = x(:);
% end
% 
% end
% %--------------------------------------------------------------------


%--------------------------------------------------------------------------
function [Ck, LabCk] = updatecenters(ierr, Q, I)
[X,Y,C] = size(I);
indC = X*Y*(0:C-1);

nk = length(ierr);
Ck = ones(nk,2);
LabCk = zeros(nk,3);
for k=1:nk
    kk = ierr(k);
    [i,j] = find(Q==kk);  % or ind2sub([X Y],Q==kk);
    % Ck(k) = sub2ind([X Y], floor(sum(i)/length(i)), floor(sum(j)/length(j)));
    Ck(k,:) = [floor(sum(i)/length(i)), floor(sum(j)/length(j))];
    ij = sub2ind([X Y],i,j);
    LabCk(k,:) = mean(I(repmat(ij,[1 C])+repmat(indC,[length(ij) 1])));
end
end
