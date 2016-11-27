%% HISTDISSIMILARITY - Dissimilarity feature derived from (co)occurrence matrices.
%
%% Description
% Calculate the dissimilarity feature defined in [HSD73] from (co)occurrence
% matrices.
%
% The local dissimilarity |D| is calculated as [ST99]: 
%
% * GLCM case - with $P_{ij}$ the joint probability of the greylevel pair
%     $(i,j)$ : $D = \sum_{i,j} P_{ij} \cdot |i-j|$
% * GLSDV case - with $P_d$ the probability of the greylevel difference  
%     $d=i-j$ : $D = \sum_d P_d \cdot |d|$
%
%% Syntax
%     D = HISTDISSIMILARITY(Pij, i, j);
%     D = HISTDISSIMILARITY(Pd, d);
%
%% References
% See |HISTCONTRAST|. 
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
% <HISTIDIFFERENCE.html |HISTIDIFFERENCE|>,
% <HISTCORRELATION.html |HISTCORRELATION|>,
% <LOCALGLCM2D.html |LOCALGLCM2D|>,
% <LOCALGLOV2D.html |LOCALGLOV2D|>,
% <LOCALGLSDV2D.html |LOCALGLSDV2D|>.

%% Function implementation
function D = histdissimilarity(pIJ,I,varargin)
D = I;
if nargin>2
    D = D - varargin{1};
end;

D = abs(D);
D = sum (pIJ(:) .* D(:));
end % end of histdissimilarity
