%% GSTSMOOTH_BASE - Base function for GSTSMOOTH.
%
%% Syntax
%      [T, gx, gy, gx2, gy2, gxy] = GSTSMOOTH_BASE(I, rho, sig, der, int, ...
%                          samp, hsize, gn, tn, thez, sigt);
%
%% Acknowledgment
% This source code uses the coherence filtering toolbox of D.Kroon.
%
%% See also
% Related:
% <GSTSMOOTH.html |GSTSMOOTH|>,
% <HESSMOOTH_BASE.html |HESSMOOTH_BASE|>.
% Called:
% <GRDSMOOTH_BASE.html |GRDSMOOTH_BASE|>,
% <GRD2GST_BASE.html |GRD2GST_BASE|>,
% <../../statistics/html/UPSCALEXY_BASE.html |UPSCALEXY_BASE|>.

%% Function implementation
function [T,gx,gy,gx2,gy2,gxy] = ...
    gstsmooth_base(I, rho, sig, der, int, samp, hsize, gn, tn, thez, sigt)

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
% * upsampling: interpolation following the recommandation of [Koth03]

if samp>1
    S =[samp samp]; % interpolation of a factor samp in X- and Y-directions
    % perform interpolation channel by channel
    I = upscalexy_base(I, S, 'linear');
end

%%
% * computation of the directional derivatives through Gaussian smoothing
% and differentiation

sig = sig * samp;

%%
% compute the gradient (gy: vertical, gx: horizontal) for each channel
[gy,gx] = grdsmooth_base(I, sig, der, hsize{1}, 'ij');
% same as computing first: [gx,gy]=grdsmooth(I,sigma,p.der,hsize,'xy');
% and then take the vector orthogonal to the gradient: tmp=gx; gx=gy; gy=-tmp;
gy = -gy;
% note that the output directional derivatives have size [X Y C]    

%% 
% * renormalization of the gradient by its norm, channel by channel

if gn
    d = sqrt( gx.^2 + gy.^2 );
    d(d<eps) = 1;
    gx = gx ./ d; gy = gy ./ d;
end

%% 
% * estimation of the gradient structure tensor

rho = rho * samp;
[gx2,gy2,gxy] = grd2gst_base(gx, gy, rho, int, hsize{2}, samp, tn, thez, sigt, []);

%% 
% * output: assemble the output tensor 

% T = zeros(X,Y,2,2);
% T(:,:,1,1) = gx2(1:samp:end,1:samp:end);
% T(:,:,2,2) = gy2(1:samp:end,1:samp:end);
% T(:,:,1,2) = gxy(1:samp:end,1:samp:end);

T = zeros(samp*X,samp*Y,2,2);
T(:,:,1,1) = gx2;
T(:,:,2,2) = gy2;
T(:,:,1,2) = gxy;

T(:,:,2,1) = T(:,:,1,2);

% if nargout>=2
%     gx = gx(1:samp:end,1:samp:end,:);
%     gy = gy(1:samp:end,1:samp:end,:);
% end

end % end of gstsmooth_base



