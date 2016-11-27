%% HISTOINTEGRALJOINT_BASE - Joint integral histogram of an image.
%
%% Description
% Compute the joint integral histograms defined in [ZLLG10].
%
%% Syntax
%     IH = HISTOINTEGRALJOINT_BASE(I, f, nbin, bound);
%
%% References
% [VJ01]  P. Viola and M. Jones, "Robust Real-Time Face Detection", 
%      International  Journal of Computer Vision, 57(2):137-154, 2004.
%
% [Pori05]  F. Porikli: "Integral histogram: a fast way to extract histogram
%      features", Proc. IEEE CVPR, pp. 829?836, 2005.
%
% [PT06]  F. Porikli and O. Tuzel: "Fast construction of covariance matrices
%      for arbitrary size image windows", Proc. IEEE ICIP, 2006.
%
% [ZLLG10]  K. Zhang, G. Lafruit, R. Lauwereins and L. Van Gool: "Joint 
%      integral histograms and its application in stereo matching", Proc. 
%      IEEE ICIP, pp. 817-820, 2010.
%
%% See also
% Related:
% <HISTOINTEGRAL1D_BASE.html |HISTOINTEGRAL1D_BASE|>,
% <HISTOINTEGRAL2D_BASE.html |HISTOINTEGRAL2D_BASE|>.
% Called:
% <INTEGRALIMAGE.html |INTEGRALIMAGE|>.

%% Function implementation
function IH = histointegraljoint_base(I, f, nbin, bound)

[X,Y] = size(I(:,:,1));       

if nargin==4,
    imax = bound(2); imin = bound(1);
else
    imax = max(I(:)); imin = min(I(:));
    if nargin<3
        nbin = round(imax - imin) + 1;
    end
end

% quantize the input image
I = floor((nbin-1) * (I - imin) / (imax - imin));

% optimized, fully vectorial code
ind = (1:X*Y)';
IH = zeros(X, Y, nbin);

IH(ind + X*Y*I(:)) = f(:);
%IH(ind + X*Y*I(:)) = 1;
%F = cumsum(f,2);
%IH = IH .* repmat(F, [1 1 nbin]);

% IH = cumsum(cumsum(IH,2),1);
IH = integralimage(IH);

end % end of histointegraljoint_base
