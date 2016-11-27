%% HISTENTROPY10 - Entropy feature derived from (co)occurrence matrices.
% 
%% Description
% Compute the entropy feature defined in [HSD73] from (co)occurrence matrices
% following the same formula as the function |HISTENTROPY2|, but using the 
% the $\log$ in base 10 instead. 
%
%% Syntax
%     E = HISTENTROPY10(Pij);
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
% <HISTMAXIMUM.html |HISTMAXIMUM|>,
% <HISTMEAN.html |HISTMEAN|>,
% <HISTDISSIMILARITY.html |HISTDISSIMILARITY|>,
% <HISTIDIFFERENCE.html |HISTIDIFFERENCE|>,
% <HISTCORRELATION.html |HISTCORRELATION|>,
% <LOCALGLCM2D.html |LOCALGLCM2D|>,
% <LOCALGLOV2D.html |LOCALGLOV2D|>,
% <LOCALGLSDV2D.html |LOCALGLSDV2D|>.

%% Function implementation
function E = histentropy10(pIJ, varargin)
pIJ(pIJ==0) = 1;
E = - sum(log10(pIJ(:)) .* pIJ(:));
end % end of histentropy10
