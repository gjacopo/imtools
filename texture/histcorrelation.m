%% HISTCORRELATION - Correlation feature derived from (co)occurrence matrices.
%
%% Description
% Calculate the correlation feature defined in [HSD73] from  (co)occurrence
% matrices.
%
% This feature measures the correlation |C| of pixel pairs on greylevels
% as [ST99]:
%
% * GLCM case - with the local statistics $\mu_i = \sum_{i,j} i \cdot P_{ij}$,
%     $\mu_j = \sum_{i,j} j \cdot P_{ij}$, $\sigma_i = \sum_{i,j} (i-\mu_i)^2 \cdot P_{ij}$ 
%     and $\sigma_j = \sum_{i,j} (j-\mu_j)^2 \cdot P_{ij}$ based on the  
%     greylevel joint probability $P_{ij}$ of the greylevel pair $(i,j)$ :
% $$ 
%     C = \sum_{i,j} ((i-\mu_i) \cdot (j-\mu_j) \cdot P_{ij}) ) / (\sigma_i \cdot \sigma_j)
%
% * GLSDV case - with the local statistics $\mu_d = \sum_d d \cdot P_d$,
%     $\mu_s = \sum_s s \cdot P_s$ and
% $$ 
%     \sigma = \sum_s (s-\mu_s)^2 \cdot P_s + \sum_d (d-\mu_d)^2 \cdot P_d
% $$ 
%     based on the probabilities $P_s$ and $P_d$ of the greylevel sum $s$ and
%     difference $d$ respectively :
% $$ 
%     C = \frac{(\sum_s(s-\mu_s)^2 \cdot P_s - \sum_d (d-\mu_d)^2 \cdot P_d)}{\sigma^2}
% $$      
%
%% Syntax
%     C = HISTCORRELATION(pD, D);
%     C = HISTCORRELATION(pIJ, I, J);
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
% <HISTIDIFFERENCE.html |HISTIDIFFERENCE|>,
% <LOCALGLCM2D.html |LOCALGLCM2D|>,
% <LOCALGLOV2D.html |LOCALGLOV2D|>,
% <LOCALGLSDV2D.html |LOCALGLSDV2D|>.

%% Function implementation
function C = histcorrelation(pIJ,I,varargin)

muI = sum(I .* pIJ);
sI = sum((I-muI).^2 .* pIJ);

if nargin>2
    J = varargin{1}; 
    muJ = sum(J .* pIJ);
    sJ = sum((J-muJ).^2 .* pIJ);
else % J = I = D
    J = I;
    muJ = muI;
    sJ = sI;
end

% Correlation = sum_i( sum_j( ((ij)p(i,j) - u_x.u_y) / (s_x.s_y) ) ) (p[2])
% C = sum (I .* J .* pIJ - muI .* muJ) / sqrt(sI * sJ);
% Correlation = sum_i( sum_j( (i - u_i)(j - u_j)p(i,j)/(s_i.s_j) ) ) (m)
C = sum ((I - muI) .* (J - muJ) .* pIJ) / sqrt(sI * sJ);

end % end of histcorrelation
