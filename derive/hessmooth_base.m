%% HESSMOOTH_BASE - Base function for HESSMOOTH.
%
%% Syntax
%     [H, gx, gy] = HESSMOOTH_BASE(I, rho, sig, der, int, samp, hsize, gn, tn, thez, sigt);
%
%% See also  
% Related:
% <HESSMOOTH.html |HESSMOOTH|>,
% <GSTSMOOTH_BASE.html |GSTSMOOTH_BASE|>.
% Called:
% <GRDSMOOTH_BASE.html |GRDSMOOTH_BASE|>,
% <GRD2HESS_BASE.html |GRD2HESS_BASE|>,
% <../../statistics/html/UPSCALEXY_BASE.html |UPSCALEXY_BASE|>.

%% Function implementation
function [H,gx,gy,gxx,gyy,gxy] = ...
    hessmooth_base(I,rho, sig, der, int, samp, hsize, gn, tn, thez, sigt)

%%
% setting internal variables

[X,Y,C] = size(I);                                                     %#ok   

%%
% define the filter window size: it is passed in order to constrain the
% smoothing and integration window size
if ~iscell(hsize)
    if length(hsize)==2,   j=2;
    else                   j=1;  end
    tmp = hsize;   hsize = cell(2,1);
    if ~isempty(tmp)
        hsize{1} = tmp(1);  hsize{2} = tmp(j);
    end
    % when hsize=[], then it is set a a cell of empty matrices
end

%%
% * computation of the directional derivatives through Gaussian smoothing and
% differentiation

%%
% 1st order derivation: compute the gradient (gy: vertical, gx: horizontal)
% for each channel
[gy,gx] = grdsmooth_base(I, sig, der, hsize{1}, 'ij'); 
gy = -gy;

%% 
% * resampling 

%%
% interpolation following the recommandation of [Koth03]
if samp>1
    S =[samp samp]; % interpolation of a factor samp in X- and Y-directions
    % perform interpolation channel by channel
    gx = upscalexy_base(gx, S, 'linear');
    gy = upscalexy_base(gy, S, 'linear');
end

%% 
% * renormalization of the gradient by its norm, channel by channel
if gn
    d = sqrt( gx.^2 + gy.^2 );
    d(d<eps) = 1;
    gx = gx ./ d; gy = gy ./ d;
end

%% 
% * estimation of the smoothed Hessian estimation
[gxx,gyy,gxy] = grd2hess_base(gx, gy, rho, der, int, ...
    hsize{2}, samp, tn, thez, sigt);

%% 
% * output: assemble the output Hessian tensor

H(:,:,1,1) = gxx(1:samp:end,1:samp:end);
H(:,:,2,2) = gyy(1:samp:end,1:samp:end);
H(:,:,1,2) = gxy(1:samp:end,1:samp:end);
H(:,:,2,1) = H(:,:,1,2);

end % end of hessmooth_base



