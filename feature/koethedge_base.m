%% KOETHEDGE_BASE - Base function for KOETHEDGE.
% 
%% Syntax
%     [emap, cmap, T] = ...
%         KOETHEDGE_BASE(I, rho, sigma, der, int, map, samp, v1, v2, radius);
%
%% See also
% Related:
% <KOETHEDGE.html |KOETHEDGE|>,
% <CANNYEDGE_BASE.html |CANNYEDGE_BASE|>,
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>,
% <CONGRUENCYEDGE_BASE.html |CONGRUENCYEDGE_BASE|>,
% <COMPASSEDGE_BASE.html |COMPASSEDGE_BASE|>,
% <ANISOEDGE_BASE.html |ANISOEDGE_BASE|>,
% <ELDERZUCKEREDGE_BASE.html |ELDERZUCKEREDGE_BASE|>,
% <ROTHWELLEDGE_BASE.html |ROTHWELLEDGE_BASE|>,
% <SDGDEDGE_BASE.html |SDGDEDGE_BASE|>,
% <PETROUEDGE_BASE.html |PETROUEDGE_BASE|>.
% Called: 
% <../../derive/html/GSTSMOOTH_BASE.html |GSTSMOOTH_BASE|>,
% <../../derive/html/GSTDECOMP.html |GSTDECOMP|>,
% <../../derive/html/GSTFEATURE_BASE.html |GSTFEATURE_BASE|>,
% <../../statistics/html/FINDLOCALEXTREMA_BASE.html |FINDLOCALEXTREMA_BASE|>,
% <CANNYEDGEMAP_BASE.html |CANNYEDGEMAP_BASE|>.

%% Function implementation
function [emap, cmap, T] = ...
    koethedge_base(I, rho, sigma, der, int, map, samp, v1, v2, radius )

% note : samp=2 and int='ani' corresponds to the method proposed by Koethe

% only parameters set to constant values
thez =8;  sigt = .4; gn = false; tn = false; hsize = [];

%%
% compute the gradient structure tensor
[T, ~, ~, gx2, gy2, gxy] = ...
 gstsmooth_base(I, rho, sigma, der, int, samp, hsize, gn, tn, thez, sigt);

%%
% perform the eigen decomposition
[l1,l2,e1,~] = gstdecomp(T);                                         
% % compute Tedge = (l1-l2) e1 . e1^T
% Tedge = zeros(size(T));
% l1 = l1 - l2;
% Tedge(:,:,1,1) = l1 .* e1(:,:,1) .^2;
% Tedge(:,:,2,2) = l1 .* e1(:,:,2) .^2;
% Tedge(:,:,1,2) = l1 .* e1(:,:,1) .* e1(:,:,2);
% Tedge(:,:,1,2) = Tedge(:,:,2,1);
%
% % compute Tjunction = l2 * I
% Tjunction = zeros(size(T));
% Tjunction(:,:,1,1) = l2;
% Tjunction(:,:,2,2) = l2;

%%
% ... and feed it to a Canny-like algorithm: non maxima suppression and
% hystheresis thresholding
switch map
    
    case 'canny'
        % transform the tensor back into a gradient-like vector
        mag = gstfeature_base(gx2, gy2, gxy, 'norm', v1, [], []);
        % compute the double angle of the gradient tensor
        or = gstfeature_base(gx2, gy2, gxy, 'orient', [], [], []);
        emap = cannyedgemap_base(e1(:,:,1), e1(:,:,2), 'matlab', mag, or, v2, []);
        
    case 'rothwell' % subpixel detection
        mag = sqrt(l1-l2);
        emap = rothwelledge_base(mag.*e1(:,:,1), mag.*e1(:,:,2), v1, v2);
        
end

%%
% detect corners
cmap = 2*l2; % trace(Tjunction)
cmap(abs(cmap)<eps) = 0;
cmap = findlocalextrema_base(cmap, 'max', samp*radius, 8, []);

end % end of koethedge_base