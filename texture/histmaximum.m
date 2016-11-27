%% HISTMAXIMUM - Maximum feature derived from (co)occurrence matrices.
%
%% Description
% Calculate the maximum feature defined in [HSD73] from (co)occurrence matrices.
%
% This feature measures the local maximum |M| in an image as [ST99]: 
%
% * GLCM case - with $P_{ij}$ the joint probability  of the greylevel pair
%     $(i,j)$ : $M = \max_i  P_{ij}$
% * GLOV case - with the probability $P_i$ of the greylevel $i$ : $M = \max_i P_i$
%
%% Syntax
%     M = HISTMAXIMUM(Pij);
%     M = HISTMAXIMUM(Pi);
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
% <HISTMEAN.html |HISTMEAN|>,
% <HISTDISSIMILARITY.html |HISTDISSIMILARITY|>,
% <HISTIDIFFERENCE.html |HISTIDIFFERENCE|>,
% <HISTCORRELATION.html |HISTCORRELATION|>,
% <LOCALGLCM2D.html |LOCALGLCM2D|>,
% <LOCALGLOV2D.html |LOCALGLOV2D|>,
% <LOCALGLSDV2D.html |LOCALGLSDV2D|>.

%% Function implementation
function M = histmaximum(pIJ, varargin)
M = max (pIJ(:));
end % end of histmaximum
