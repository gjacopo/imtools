%% HISTMEAN - Mean feature derived from (co)occurrence matrices.
%
%% Description
% Calculate the mean feature defined in [HSD73] from (co)occurrence matrices.
% 
% This feature measures the average |A| of the greylevel within an image 
% as [ST99]:
%
% * GLCM case - with $P_{ij}$ the greylevel joint probability of the greylevel 
%     pair $(i,j)$ : $A = \sum_{i,j} P_{ij}$ 
%
% * GLOV case - with $P_i$ the probability of the greylevel $i$ : $A = \sum_i P_i$ 
%
%% Syntax
%     M = HISTMEAN(Pij);
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
% <HISTDISSIMILARITY.html |HISTDISSIMILARITY|>,
% <HISTIDIFFERENCE.html |HISTIDIFFERENCE|>,
% <HISTCORRELATION.html |HISTCORRELATION|>,
% <LOCALGLCM2D.html |LOCALGLCM2D|>,
% <LOCALGLOV2D.html |LOCALGLOV2D|>,
% <LOCALGLSDV2D.html |LOCALGLSDV2D|>.

%% Function implementation
function M = histmean(pIJ,varargin)
M = mean (pIJ(:));
end % end of histmean
