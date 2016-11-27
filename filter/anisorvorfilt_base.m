%% ANISORVORFILT_BASE - Base function for ANISORVORFILT. NOT FINISHED IMPLEMENTING. 
%
%% Description
% Morphological anisotropic filtering of connected structures using the
% orientation information derived from the anisotropic continous scale model
% of [BBW07] (numerical scheme therein) and within the Voronoi diagram of
% the connected structures. 
%
%% Syntax
%     u = ANISORVORFILT_BASE(I, delta, nu, rho, sigma, der, int, samp);
%
%% References
% [BBW07]  M. Breuﬂ, B. Burgeth and J. Weickert: "Anisotropic continuous-
%      scale morphology", Proc. IbPRIA, LNCS 4478, pp. 515-522, Springer, 
%      2007.	
%      <http://www.springerlink.com/content/1hm264w86111m148/>
%
% [PP08]  G. Papari and N. Petkov: "Adaptive pseudo dilation for Gestalt 
%      edge grouping and contour detection", IEEE Trans. on Image Processing,
%      17(10):1950-1962, 2008.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4623243&tag=1>
%
%% See also
% Related:
% <ANISOR.html |ANISOR|>,
% <ANISORFILT.html |ANISORFILT|>,
% <../../derive/html/GSTSMOOTH.html |GSTSMOOTH|>,
% <matlab:webpub(whichpath('VORONOI')) |VORONOI|>.
% Called: 
% <ANISOR_BASE.html |ANISOR_BASE|>,
% <../../propagation/html/POTENTIAL2FRONT.html |POTENTIAL2FRONT|>,
% <matlab:webpub(whichpath('BWCONNCOMP')) |BWCONNCOMP|>,
% <matlab:webpub(whichpath('LABELMATRIX')) |LABELMATRIX|>,
% <matlab:webpub(whichpath('IMDILATE')) |IMDILATE|>,
% <matlab:webpub(whichpath('IMERODE')) |IMERODE|>.

%% Function implementation
function u = anisorvorfilt_base(I, delta, nu, rho, sigma, der, int, samp) %#ok  

error('anisorvorfilt_base:exit','not implemented yet');

%%
% setting internal variables

[X,Y] = size(I);                                                       %#ok                                            

% predefine appropriate indices
zer0 = false(X+2,Y+2);

iminusj = zer0; iminusj(1:end-2,2:end-1) = true; iminusj = iminusj(:);
iplusj = zer0; iplusj(3:end,2:end-1) = true; iplusj = iplusj(:);
ijminus = zer0; ijminus(2:end-1,1:end-2) = true; ijminus = ijminus(:);
ijplus = zer0; ijplus(2:end-1,3:end) = true; ijplus = ijplus(:);

zer0 = zeros(X*Y,1);

%%
% create labelled image

if islogical(I)
    CC = bwconncomp(I);
    L = labelmatrix(CC);
    nCC = CC.NumObjects;
else
    % see
    % http://www.mathworks.com/matlabcentral/fileexchange/26946-label-connected-components-in-2-d-array
end
  
%% 
% compute the Voronoi diagram of the labelled structures

W = ones(X, Y);

[i,j] = ind2sub([X,Y], find(L~=0));
[D, VV] = potential2front(W,  [i'; j']);
V = zeros(X,Y);
for ip = 1:length(i)
    V(VV==ip) = L(i(ip),j(ip));
end

% D = zeros(X, Y, nCC);
% for l=1:nCC
%     [i,j] = ind2sub([X,Y], find(L==l));      
%     D(:,:,l) = potential2front(W,  [i'; j']);
% end
% [dum, V] = sort(D, 3, 'ascend');
% V = V(:,:,1); % Voronoi diagram of the connected components
% figure, imagesc(V)
% V(L~=0) = L(L~=0);

% frontiers
se = strel('square',3);
dV = imdilate(V,se) - imerode(V,se);
dV = dV > 0;

%%
% perform Rouy-Turin numerical scheme like in Eq.(11)

u = zeros(size(I));
uu = double(I);

%while any(~isout)
for i=1:30
    kappa = anisor_base(uu, nu, rho, sigma, der, int, samp);
    % kappa=1;
    
    %    for l=1:nCC
    %        if isout(l),  continue;  end
    U = padarray( uu, [1 1], 'replicate', 'both');
    
    % Rouy-Tourin scheme a in Eq.(11): naive implementation
    h = 1;
    uu = uu(:) + delta * kappa(:) .* ...
        sqrt((max([zer0 (U(iplusj) - uu(:))/h (U(iminusj) - uu(:))/h],[],2)).^2 + ...
        (max([zer0 (U(ijplus) - uu(:))/h (U(ijminus) - uu(:))/h],[],2)).^2);
    % note: the sign of delta defines if the morphological operation is
    % a dilation ('+' sign) or an erosion ('-' sign)
    
    % http://cermics.enpc.fr/~forcadel/Publi/FGL.pdf
    % http://etd.lsu.edu/docs/available/etd-09152008-143521/unrestricted/Tu
    % gurlandiss.pdf
    % http://ctr.stanford.edu/ResBriefs03/herrmann1.pdf
    % http://www.mia.uni-saarland.de/Publications/pizarro-ismm09.pdf

end

end % end of anisorvorfilt_base
