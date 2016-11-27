%% ANISOR_BASE - Base function for ANISOR.
%
%% Syntax
%     kappa = ANISOR_BASE(I, mu, rho, sigma, der, int, samp);
%
%% See also
% Related:
% <ANISOR.html |ANISOR|>.
% Called:
% <../../derive/html/GSTSMOOTH_BASE.html |GSTSMOOTH_BASE|>,
% <../../derive/html/GSTDECOMP.html |GSTDECOMP|>,
% <matlab:webpub(whichpath('DOT')) |DOT|>.

%% Function implementation
function kappa = anisor_base(I, mu, rho, sigma, der, int, samp)

C = size(I,3);
% possibly multichannel when C>1...or not!!!                       
if(C>1) 
    error('anisor_base:inputerror',...
        ['TODO: method not implemented yet for multichannel images - need ' ...
        'to define an appropriate multichannel gradient for such images - a ' ...
        'possible solution can be given by Scheunders approach as we are '... 
        'interested in the orientation only (see subfunction GSTREORIENT in ' ...
        'GSTFEATURE_BASE)']);
end

%%
% compute the gradient (given by the directional derivatives) and the 
% gradient structure tensor

%%
% compute the structure tensor |S| of Eq.(8)
[S, gx, gy] = ...
    gstsmooth_base(I, rho, sigma, der, int, samp, [], false, false, 8, .4);
%  |rho|: integration scale
%  |sigma|: differentiation scale
%  |der, int|: methods used for derivating and smoothing resp. see |GRDSMOOTH|,
%   |GSTSMOOTH| and |SMOOTHFILT| helps.

%%
% note: the components of |S| and the directional derivatives |(gx,gy)| were
% computed using the same |sigma|

%% 
% perform the eigen decomposition
[~, ~, ~, v2] = gstdecomp(S); 
% in [l1, l2, v1, v2] = gstdecomp(S), we have: 
%  l1, l2 : eigenvalues of the tensor, l1>=l2. 
%  v1, v2 : corresponding eigenvectors; "v1 describes the orientation of 
%   highest contrast variation", ie. it gives the direction of maximal 
%   local change  within the window given the integration scale rho, and v2
%   is orthogonal to v1.

%%
% in the case of scalar image, l1 \approx |grad|, v1 \approx orient(grad) 
% note that however there is an undetermination in the direction (sign) of
% the eigenvectors. But it is something dealt with |eta| in the following
% (see remark below)

% figure, quiver(gx,gy),title('gx gy');
% figure, quiver(v1(:,:,1),v1(:,:,2)),title('v2');
% figure, quiver(v2(:,:,1),v2(:,:,2)),title('v2');


%% 
% compute the orientation information
nabla = cat(3, gx, gy); % this is the gradient field \nabla u_{\sigma}

%%
% compute the |eta| of Eq.(9)
eta = abs(dot(nabla, v2, 3) ./ (sqrt(dot(nabla,nabla,3) .* dot(v2,v2,3))));
eta(isnan(eta))=0; % consider the case 0/0

%%
% note: when |A| and |B| are two vector fields defined over the domain |[X,Y]| 
% (therefore, |A| - ibid with |B| - is a matrix of size |(X,Y,2)| whose X- 
% and Y-coordinates are defined in |A(:,:,1)| and |A(:,:,2)| resp. - ibid with
% |B|), |dot(A,B,3)| is the scalar product of |A| and |B| calculated all over
% the the domain |[X,Y]|. Consequently, |dot(A,A,2)| is the (squared) norm 
% of the vector field |A| defined on the domain |[X,Y]|

%%
% remark: at that point, |eta| depends only on the respective orientations of
% the vectors nabla and |v2|, not on their directions

%%
% compute the kappa of Eq.(10), where |mu| is a non-negative integer which 
% influences the propagation of the structures in |v2| direction
kappa = eta .^ mu;

end % end of anisor_base
