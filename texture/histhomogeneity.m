%% HISTHOMOGENEITY - Homogeneity feature derived from (co)occurrence matrices.
%
%% Description
% Calculate the homogeneity feature defined in [HSD73] from (co)occurrence 
% matrices.
%
% This feature measures the local homogeneity |H| in an image as [ST99]:
%
% * GLCM case - with $P_{ij}$ the joint probability of the greylevel pair
%     $(i,j)$ :  
% $$ 
%     H = \sum_{i,j}  \frac{P_{ij}}{(1+(i-j)^2)}
% $$ 
%
% * GLSDV case - with $P_d$ the probability of the greylevel difference  
%     $d=i-j$ :
% $$ 
%     H = \sum_d \frac{P_d}{(1+d^2)}
% $$ 
%
%% Syntax
%     H = HISTHOMOGENEITY(Pij, i, j);
%     H = HISTHOMOGENEITY(Pd, d);
%
%% References
% See |HISTCONTRAST|. 
%
%% See also
% <HISTCONTRAST.html |HISTCONTRAST|>,
% <HISTENERGY.html |HISTENERGY|>,
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
function H = histhomogeneity(pIJ,I,varargin)

H = I;
if nargin>2
    H = H - varargin{1};
end;

H = (1 + H .* H);
H = sum (pIJ(:) ./ H(:));
end % end of histhomogeneity
