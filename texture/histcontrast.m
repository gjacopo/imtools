%% HISTCONTRAST - Contrast feature derived from (co)occurrence matrices.
%
%% Description
% Calculate the contrast feature defined in [HSD73] from (co)occurrence
% matrices.
%
% This feature measures the local contrast |C| in an image as [ST99]: 
%
% * GLCM case - with $P_{ij}$ the joint probability of the greylevel pair
%     $(i,j)$ :
% $$ 
%     C = \sum_{i,j} (i-j)^2 \cdot P_{ij}
% $$ 
%
% * GLSDV case - with $P_d$ the probability of the greylevel difference 
%     $d=i-j$ :
% $$ 
%     C = \sum_d d^2 \cdot P_d
% $$ 
%
%% Syntax
%     C = HISTCONTRAST(Pij, i, j);
%     C = HISTCONTRAST(Pd, d);
%
%% References
% [HSD73] R.M. Haralick, K. Shanmugam, and I. Dinstein: "Textural features
%      for image classification", IEEE Trans. Systems, Man and Cybernetics,
%      3(6):610-621, 1973.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4309314>
%
% [Hara79] R.M. Haralick: "Statistical and structural approaches to texture",
%      Proceedings of IEEE, 67:786-804, 1979. 
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1455597>
%
% [ST99] L. Soh and C. Tsatsoulis: "Texture analysis of SAR sea ice imagery
%      using gray level co-occurrence matrices", IEEE Trans. on Geoscience
%      and Remote Sensing, 37(2), 1999.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=752194>
% 
%% See also
% Related:
% <HISTENERGY.html |HISTENERGY|>,
% <HISTHOMOGENEITY.html |HISTHOMOGENEITY|>,
% <HISTVARIANCE.html |HISTVARIANCE|>,
% <HISTENTROPY2.html |HISTENTROPY2|>,
% <HISTENTROPY10.html |HISTENTROPY10|>,
% <HISTMAXIMUM.html |HISTMAXIMUM|>,
% <HISTMEAN.html |HISTMEAN|>,
% <HISTDISSIMILARITY.html |HISTDISSIMILARITY|>,
% <HISTIDIFFERENCE.html |HISTIDIFFERENCE|>,
% <HISTCORRELATION.html |HISTCORRELATION|>,
% <LOCALGLCM2D.html |LOCALGLCM2D|>,
% <LOCALGLOV2D.html |LOCALGLOV2D|>,
% <LOCALGLSDV2D.html |LOCALGLSDV2D|>.

%% Function implementation
function C = histcontrast(pIJ,I,varargin)

C = I;
if nargin>2
    C = C - varargin{1};
end;

C = C .* C;
C = sum (C(:) .* pIJ(:));
end % end of histcontrast
