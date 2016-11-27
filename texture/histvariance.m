%% HISTVARIANCE - Variance feature derived from (co)occurrence matrices.
%
%% Description
% Calculate the variance feature defined in [HSD73] from (co)occurrence 
% matrices.
% 
% This feature measures the variation |V| of greylevel distribution as [ST99]:
%
% * GLCM case - with the local statistics $\mu_i = \sum_{i,j} i \cdot P_{ij}$
%     and $\mu_j = \sum_{i,j} j \cdot  P_{ij}$, based on the greylevel joint 
%     probability $P_{ij}$ of the greylevel pair $(i,j)$ :
% $$ 
%     V = \sum_{i,j}(i \cdot j \cdot P_{ij} - \mu_i \cdot \mu_j)
% $$ 
%
% * GLSDV case - with $\mu_d = \sum_d d \cdot P_d$ the local mean based on 
%     the probability $P_d$ of the greylevel difference $d=i-j$ :  
% $$
%     V = \sum_d d^2 \cdot P_d - \mu_d^2
% $$ 
%
%% Syntax
%     V = HISTVARIANCE(Pij, i, j);
%     V = HISTVARIANCE(Pd, d);
%
%% References
% See |HISTCONTRAST|. 
% <mailto:grazzja@lanl.gov J.Grazzini> (ISR-2/LANL)
%
%% See also
% Related:
% <HISTCONTRAST.html |HISTCONTRAST|>,
% <HISTENERGY.html |HISTENERGY|>,
% <HISTHOMOGENEITY.html |HISTHOMOGENEITY|>,
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
function V = histvariance(pIJ, I, varargin)

muI = sum( I(:) .* pIJ(:)); 

if nargin>2
    J = varargin{1};
    muJ = sum( J(:) .* pIJ(:));
    V = sum (I(:) .* J(:) .* pIJ(:) - muI .* muJ);
    
else % J = I = D
    V = sum (I(:).^2 .* pIJ(:) - muI.^2);
end

end % end of histvariance
