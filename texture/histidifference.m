%% HISTIDIFFERENCE - Inverse difference feature derived from (co)occurrence matrices.
%
%% Description
% Calculate the inverse difference feature defined in [HSD73] from (co)occurrence 
% matrices.
%
% This feature measures the local inverse difference |D| in an image as [ST99]: 
%
% * GLCM case - with $P_{ij}$ the joint probability of the greylevel pair 
%     $(i,j)$ :  
% $$ 
%     D = \sum_{i,j} \frac{P_{ij}}{(1+|i-j|)}
% $$ 
%
% * GLSDV case - with $P_d$ the probability of the greylevel difference
%     $d=i-j$ :
% $$ 
%     D = \sum_d \frac{P_d}{(1+|d|)}
% $$ 
% 
%% Syntax
%     D = HISTIDIFFERENCE(Pd, d);
%     D = HISTIDIFFERENCE(Pij, i, j);
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
% <HISTDISSIMILARITY.html |HISTDISSIMILARITY|>,
% <HISTCORRELATION.html |HISTCORRELATION|>,
% <LOCALGLCM2D.html |LOCALGLCM2D|>,
% <LOCALGLOV2D.html |LOCALGLOV2D|>,
% <LOCALGLSDV2D.html |LOCALGLSDV2D|>.

%% Function implementation
function D = histidifference(pIJ,I,varargin)
D = I;
if nargin>2
    D = D - varargin{1};
end;

D = 1 + abs(D);
D = sum (pIJ(:) ./ D(:));

end % end of histidifference
