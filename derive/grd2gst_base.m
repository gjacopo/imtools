%% GRD2GST_BASE - Base function for GRD2GST.
% 
%% Syntax
%     [gx2, gy2, gxy] = GRD2GST_BASE(gx, gy, rho, int, hsize, ...
%                                    samp, tn, thez, sigt, eign);
%     [gx2, gy2, gxy, mag, or] = GRD2GST_BASE(gx, gy, rho, int, hsize, ...
%                                           samp, tn, thez, sigt, eign);
%
%% See also  
% Related:
% <GRD2GST.html |GRD2GST|>.
% Called:
% <GSTDECOMP.html |GSTDECOMP|>,
% <GSTFEATURE_BASE.html |GSTFEATURE_BASE|>,
% <../../filter/html/SMOOTHFILT_BASE.html |SMOOTHFILT_BASE|>.

%% Function implementation
function [gx2, gy2, gxy, varargout] = ...
    grd2gst_base(gx, gy, rho, int, hsize, samp, tn, thez, sigt, eign)

[X Y C] = size(gx);                                                    %#ok       

%%
% * compute unsmoothed tensor

%%
% compute the component of J(grad u_sigma) while taking into account the
% multichannel information by averaging over the different channels
gx2 = mean(gx.^2, 3);
gxy = mean(gx.*gy, 3);
gy2 = mean(gy.^2, 3);
%%
% following the approach of [Scheun03], average the gradient over the
% channels, which is nothing else than the gradient of the average of the
% channels

%%
% * renormalization of the tensor

if tn
    sigma0 = 16; % initial isotropy
    eta0 = 12; % intial anisotropy
    [l1,l2,e1,e2] = gstdecomp(gx2, gy2, gxy); 
    l1 = l1*0 + sqrt(sigma0*eta0);
    l2 = l2*0 + sqrt(sigma0/eta0);
    [gx2, gy2, gxy] = gstdecomp(l1,l2,e1,e2);
end

%%
% * linear or nonlinear spatial integration

%%
% compute the direction of smoothing (along edges) in the anisotropic case
if ischar(int) && any(strcmp(int,{'ani','hourglass'}))
    theta = mod(atan2(mean(gy,3),mean(gx,3)),pi);
else
    theta =[];
end

%%
% perform tensor smoothing: isotropic (classical) or anisotropic (improved)
% spatial averaging of the tensor and of the directional derivatives

if rho~=0 && ~(islogical(int) && ~int)
    gx2 = smoothfilt_base( gx2, rho, int, hsize, samp, thez, sigt, theta );
    gy2 = smoothfilt_base( gy2, rho, int, hsize, samp, thez, sigt, theta );
    gxy = smoothfilt_base( gxy, rho, int, hsize, samp, thez, sigt, theta );
end

%%
% * output
if nargout >= 4
    % norm based on the eigenvalues
    varargout{1} = gstfeature_base(gx2, gy2, gxy, 'norm', eign, [], []);
    %  double angle of the gradient tensor
    if nargout==5
        varargout{2} = gstfeature_base(gx2, gy2, gxy, 'orient', [], [], []);
    end
end

end % end of grd2gst_base
