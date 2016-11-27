%% GRD2HESS_BASE - Base function for GRD2HESS.
% 
%% Syntax
%        [gxx, gyy, gxy] = ...
%    GRD2HESS_BASE(gx, gy, rho, int, hsize, samp, tn, thez, sigt);
%
%% See also  
% Related:
% <GRD2HESS.html |GRD2HESS|>.
% Called:
% <GRDSMOOTH_BASE.html |GRDSMOOTH_BASE|>,
% <../../filter/html/SMOOTHFILT_BASE.html |SMOOTHFILT_BASE|>.

%% Function implementation
function [gxx, gyy, gxy] = ...
    grd2hess_base(gx, gy, rho, der, int, hsize, samp, tn, thez, sigt)

%%
% * compute (unsmoothed) 2nd order derivative
[gxy,gxx] = grdsmooth_base(gx, 0, der, hsize, 'ij');
gxy = -gxy;
[gyy,gyx] = grdsmooth_base(gy, 0, der, hsize, 'ij');
gyy = -gyy;
gxy = (gxy+gyx) / 2;

%% 
% * combine the vectorial components (non standard operation for Hessian!)
% sum over the channels using di Zenzo approach
gxx = mean(gxx, 3);
gxy = mean(gxy, 3);
gyy = mean(gyy, 3);

%% 
% * renormalization of the Hessian tensor
if tn
    nn = gx.^2 + gy.^2;
    nn = sqrt(eps+nn);
    gxx = gxx ./ nn;
    gxy = gxy ./ nn;
    gyy = gyy ./ nn;
end

%% 
% * linear or nolinear spatial integration

%%
% compute the direction of smoothing (along edges) in the anisotropic case
if any(strcmp(int,'ani'))
    theta = mod(atan2(mean(gy,3),mean(gx,3)),pi);
else
    theta =[];
end

%% 
% perform tensor smoothing: isotropic (classical) or anisotropic (improved)
% spatial averaging of the tensor and of the directional derivatives
if rho~=0 && ~(islogical(int) && ~int)
    gxx = smoothfilt_base( gxx, rho, int, hsize, samp, thez, sigt, theta );
    gyy = smoothfilt_base( gyy, rho, int, hsize, samp, thez, sigt, theta );
    gxy = smoothfilt_base( gxy, rho, int, hsize, samp, thez, sigt, theta );
end
    
end % end of grd2hess_base

