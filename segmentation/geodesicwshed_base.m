%% GEODESICWHSED_BASE - Base function for GEODESICWSHED.
%
%% Syntax
%   [D, Q, markers] = geodesicwshed_base(I, M, method, ...
%                                  rho, sigma, a, der, int, samp, eign);
%
%% See also
% Related:
% <GEODESICWSHED.html |GEODESICWSHED|>.
% Called:
% <../../propagation/html/IM2POTENTIAL.html |IM2POTENTIAL|>,
% <../../propagation/html/POTENTIAL2FRONT.html |POTENTIAL2FRONT|>,
% <../../statistics/html/FINDLOCALEXTREMA_BASE.html |FINDLOCALEXTREMA_BASE|>,
% <../../derive/html/GSTDECOMP.html |GSTDECOMP|>.

%% Function implementation
function [DD, QQ, M] = ...
    geodesicwshed_base(I, M, method, rho, sigma, a, der, int, samp, eign)

%%
% setting internal variables

if isscalar(M)
    % the window of analysis: should be odd
    M = round( (M-1)/2 )*2 + 1;
    % center of the window
    pad = (M-1)/2;
    % pad = max([(winsize-1)/2 ceil(3*sigma) ceil(3*rho)])
    
    % padd the original image for convenient processing
    I = padarray(I,[pad pad],'replicate','both');

else
    pad = 0;
end

%%
% compute once for all the potential function based on the gradient and/or
% the gradient structure tensor of the multispectral input image and define 
% the propagation function: it is either a scalar field (isotropic case) 
% providing with the speed (strenght) for the propagation (the higher the 
% value at one point, the faster the propagation through this point) or a 
% tensor field (anisotropic case) providing with a speed and a direction for
% the propagation.
       
[T, L] = im2potential(I, method, a, rho, sigma, der, int, samp, eign);
[~, l2, e1, e2] = gstdecomp(T);

%%
% possibly define the set of markers as the minima of the 'gradient' image

if isscalar(M)
    M = findlocalextrema_base(L, 'min', M);
    M = bwulterode(M);  % note: output of BWULTERODE is always a logical array
end

if islogical(M)
    [i,j] = find(M);
    % markers must be of size 2 x nb_start_points.
    M  = [i'; j'];
end


%% 
% perform geodesic propagation from the markers

DD = Inf(size(I));
QQ = ones(size(I));

for m=1:size(M,2)
    A = distpt(I,I(M(1,m),M(2,m)));
    l1 = 1 ./ (1+A);
    l2 = l2 ./ (1+A);
    TT = gstdecomp(l1,l2,e1,e2);
    %TT =  1 ./ (eps+A);
    %min(TT(:))
    %max(TT(:))
    D = potential2front(TT, M(:,m));
    [DD, Q] = min(cat(3,D,DD), [], 3);
    QQ(Q==1) = m;
end
    
%   %----------------------------------------------------------------------
    function D = distpt(I,I0)
        [X,Y,C] = size(I);
        %D = sum((I - repmat(I0,[X,Y,C])).^2,3);
        D = sum(abs(I - repmat(I0,[X,Y,C])),3);
    end
%   %----------------------------------------------------------------------

% [DD, QQ] = potential2front(T, markers);

if pad
    DD = DD(pad+1:end-pad,pad+1:end-pad);
    QQ = QQ(pad+1:end-pad,pad+1:end-pad);
end

end % end of geodesicwshed_base

