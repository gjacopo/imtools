%% HISTENTROPY2 - Entropy feature derived from (co)occurrence matrices.
%
%% Description
% Compute the entropy feature in bits ($\log_2$ used), as defined in [HSD73]
% from (co)occurrence matrices.
%
% This feature measures the randomness |E| of greylevel distributions as [ST99]:
%
% * GLCM case - $E = - \sum_{i,j} P_{ij} \log P_{ij}$,
% * GLSDV case - $E = - \sum_d P_d \cdot \log P_d$,
% * GLOV case - $E = - \sum_i P_i \cdot \log P_i$,
%
% where the $\log$ is taken in base 2 and where the matrices (or vectors in
% the last two cases) $P_{ij}$, $P_d$ and $P_i$ represent respectively the
% probabilities of cooccurring greylevel pair $(i,j)$, of occurring greylevel 
% difference $d=i-j$ and of occurring greylevel $i$. 
%
%
%% Syntax
%     E = HISTENTROPY2(Pij);
%     E = HISTENTROPY2(Pd);
%     E = HISTENTROPY2(Pi);
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
function E = histentropy2(pIJ, varargin)
pIJ(pIJ==0) = 1; % pIJ = pIJ + (pIJ==0);
E = - sum(log2(pIJ(:)) .* pIJ(:));
end % end of histentropy2
