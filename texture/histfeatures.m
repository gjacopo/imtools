%% HISTFEATURES - Textural features derived from (co)occurrence matrices.
%
%% Description
% Calculate all together the various textural features defined in [HSD73,Hara79]
% from (co)occurrence matrices.
% See formulas in [ST99].
% 
%% Syntax
%     F = HISTFEATURES(Pij, i, j);
%     F = HISTFEATURES(Pd, d);
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
% <HISTCONTRAST.html |HISTCONTRAST|>,
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
function features = histfeatures(pIJ, I, varargin)
probIJ = pIJ(:);
diffI = I(:);
muI = sum( diffI .* probIJ);
sI = sum((I-muI).^2 .* probIJ);

if nargin>2
    J = varargin{1}; J = J(:);
    diffI = diffI - J;
    muJ = sum(J .* probIJ);
    sJ = sum((J-muJ).^2 .* probIJ);
else % J = I = D
    J = I;
    muJ = muI;
    sJ = sI;
end

% common computations

Con = diffI .* diffI;
Hom = (1 + Con);

Dis = abs(diffI);
Idiff = 1 + Dis;

% maximum
Max = max (probIJ);

% dissimilarity
Dis = sum (probIJ .* Dis);

% contrast
Con = sum (Con .* probIJ);

% inverse difference 
Idiff = sum (probIJ ./ Idiff);

% homogeneity
Hom = sum (probIJ ./ Hom);

% energy
Ene = sum(probIJ .* probIJ);

if nargin>2
% correlation
Var = sum ((I - muI) .* (J - muJ) .* probIJ) / sqrt(sI * sJ); 
% Var / sqrt(sI * sJ);
else
% variance
Var = sum (I .* J .* probIJ - muI .* muJ);
end

% entropy
probIJ(probIJ==0) = 1;
%Ent10 = - sum(log10(probIJ) .* probIJ);
Ent2 = - sum(log2(probIJ) .* probIJ);

features = [Con Ene Ent2 Hom Var Dis Idiff Max];
end % end of histfeatures

