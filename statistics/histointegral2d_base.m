%% HISTOINTEGRAL2D_BASE - 2d integral histograms of a pair of images.
% 
%% Description
% Compute the 2d integral histograms of the joint distribution of a pair of
% images.
%
%% Syntax
%     IH = HISTOINTEGRAL2D_BASE(I1, I2, nbin, bound);
%
%% References
%   [VJ01]  P. Viola and M. Jones, "Robust Real-Time Face Detection", 
%      International  Journal of Computer Vision, 57(2):137-154, 2004.
%
%   [Pori05]  F. Porikli: "Integral histogram: a fast way to extract histogram
%      features", Proc. IEEE CVPR, pp. 829?836, 2005.
%
%   [PT06]  F. Porikli and O. Tuzel: "Fast construction of covariance matrices
%      for arbitrary size image windows", Proc. IEEE ICIP, 2006.
%
%   [ZLLG10]  K. Zhang, G. Lafruit, R. Lauwereins and L. Van Gool: "Joint 
%      integral histograms and its application in stereo matching", Proc. 
%      IEEE ICIP, pp. 817-820, 2010.
%
%% Remark
% Not to be confused with the so-called "joint integral histograms" defined
% in [ZLLG10] (see also function |HISTOINTEGRALJOINT_BASE|).
%
%% See also
% Related:
% <HISTOINTEGRAL1D_BASE.html |HISTOINTEGRAL1D_BASE|>,
% <HISTOINTEGRALJOINT_BASE.html |HISTOINTEGRALJOINT_BASE|>.
% Called:
% <INTEGRALIMAGE.html |INTEGRALIMAGE|>.

%% Function implementation
function IH = histointegral2d_base(I1, I2, nbin, bound)

[X,Y] = size(I1(:,:,1));                                                   

if nargin==4,
    imax = bound(2); imin = bound(1);
else
    imax = max([I1(:); I2(:)]); imin = min([I1(:); I2(:)]);
    if nargin<3
        nbin = round(imax - imin) + 1;
    end
end

% quantize the input images
I1 = floor((nbin-1) * (I1 - imin) / (imax - imin));
I2 = floor((nbin-1) * (I2 - imin) / (imax - imin));

% compute the integral histogram
ind = (1:X*Y)';
IH = zeros(X, Y, nbin^2);
IH(ind + X*Y*(nbin*I1(:)+I2(:))) = 1;
IH = integralimage(IH);

end % end of histointegral2d_base
