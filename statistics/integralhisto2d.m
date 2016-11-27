%% INTEGRALHISTO2D - 2d integral histograms of a pair of images.
% 
%% Description
% Compute the 2d integral histograms of the joint distribution of a pair of
% images.
%
%% Syntax
%     IH = INTEGRALHISTO2D(I1, I2, nbin, bound);
%
%% References
%   [VJ01]  P. Viola and M. Jones, "Robust Real-Time Face Detection", 
%      International  Journal of Computer Vision, 57(2):137-154, 2004.
%      <http://research.microsoft.com/~viola/Pubs/Detect/violaJones_IJCV.pdf>
%
%   [Pori05]  F. Porikli: "Integral histogram: a fast way to extract histogram
%      features", Proc. IEEE CVPR, pp. 829?836, 2005.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1467353>
%
%   [PT06]  F. Porikli and O. Tuzel: "Fast construction of covariance matrices
%      for arbitrary size image windows", Proc. IEEE ICIP, 2006.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4106846>
%
%   [ZLLG10]  K. Zhang, G. Lafruit, R. Lauwereins and L. Van Gool: "Joint 
%      integral histograms and its application in stereo matching", Proc. 
%      IEEE ICIP, pp. 817-820, 2010.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5653410>
%
%% Remark
% Not to be confused with the so-called "joint integral histograms" defined
% in [ZLLG10] (see also function |HISTOINTEGRALJOINT_BASE|).
%
%% See also
% Related:
% <INTEGRALHISTO1D.html |INTEGRALHISTO1D|>,
% <INTEGRALHISTOJOINT.html |INTEGRALHISTOJOINT|>.
% Called:
% <INTEGRALIMAGE.html |INTEGRALIMAGE|>.

%% Function implementation
function IH = integralhisto2d(I1, I2, nbin, bound)

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

end % end of integralhisto2d
