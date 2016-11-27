%% INFORGBEDGE_BASE - Base function for INFORGBEDGE
%    
%% Syntax
%     [map, mag, or, rmap, rmag, ror] = ...
%               INFORGBEDGE_BASE(I, R, rho, sig, der, int);
% 
%% See also
% Related:
% <INFORGBEDGE.html |INFORGBEDGE|>,
% <KOETHEDGE_BASE.html |KOETHEDGE_BASE|>,
% <CANNYEDGE_BASE.html |CANNYEDGE_BASE|>.
% Called: 
% <../../derive/html/GSTSMOOTH_BASE.html |GSTSMOOTH_BASE|>,
% <../../derive/html/GRDSMOOTH_BASE.html |GRDSMOOTH_BASE|>,
% <CANNYEDGEMAP_BASE.html |CANNYEDGEMAP_BASE|>.

%% Function implementation
function [map, mag, or, rmap, rmag, ror] = inforgbedge_base(I, R, rho, sig, der, int)

%%
% deal with multichannel image

[T,~,~,~,gy2,gxy] = ...
    gstsmooth_base(I, rho, sig, der, int, 1, [], false, false, [], []);                                     

%%
% compute the eigen decomposition
[l1,l2,e1,e2] = gstdecomp(T);                                          %#ok                       
% D = sqrt(abs(gx2.^2 - 2*gx2.*gy2 + gy2.^2 + 4*gxy.^2));
% % compute edge orientation from first eigenvalue
% l1 = (gx2 + gy2 + D) / 2;
% % the 2nd eigenvalue would be:  l2 = (gx2 + gy2 - D) / 2;

%%
% compute the tensor norm
mag = sqrt(l1); % Koethe suggests sqrt(l1-l2);
mag = mag /max(mag(:));
% compute edge orientation (from eigenvector tangent)
or = atan2(-gxy, l1 - gy2);

%%
% feed it to a Canny-like algorithm
map = cannyedgemap_base(e1(:,:,1), e1(:,:,2), 'matlab', mag, or, [], []);

%%
% deal with reduced image

%%
% compute the single channel image gradient
[gx, gy] = grdsmooth_base(R, sig, der, [], 'ij');          

%%
% feed it to a Canny-like algorithm
[rmap, rmag, ror] = cannyedgemap_base(gx, gy, 'matlab', [], [], [], []);

end % end of inforgbedge_base