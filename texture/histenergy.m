%% HISTENERGY - Energy feature derived from (co)occurrence matrices.
%
%% Description
% Calculate the energy (2nd angular moment) feature derived from (co)occurrence 
% matrices as defined in [HSD73].
%
% This feature measures the occurrence |E| of repeated pairs within an image
% as [ST99]:
%
% * GLCM case - with $P_{ij}$ the joint probability of the greylevel pair 
%     $(i,j)$ : $E = \sum_{i,j} (P_{ij})^2$ 
%
% * GLSDV case - with $P_d$ the probability of the greylevel difference 
%     $d=i-j$ : $E = \sum_d (P_d)^2$ 
%
% * GLOV case - with $P_i$ the probability of the greylevel $i$:
%     $E = \sum_i (P_i)^2$ 
%
%% Syntax
%     E = HISTENERGY(Pij);
%     E = HISTENERGY(Pd);
%     E = HISTENERGY(Pi);
%
%% References
% See |HISTCONTRAST|. 
%
%% See also
% Related:
% <HISTCONTRAST.html |HISTCONTRAST|>,
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
function E = histenergy(pIJ,varargin)
% if nargin>1, nothing end
E = sum(pIJ(:) .* pIJ(:));
end % end of histenergy
