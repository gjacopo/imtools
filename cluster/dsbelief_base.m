function dsbelief_base()
%
% 
%      [ZTB08] Zlatoff, N. and Tellez, B. and Baskurt, A. (2008): "Combining 
%         local belief from low-level primitives for perceptual grouping", 
%       Pattern Recognition 41:1215-1229,

% similarity

% closure/compactness

end


function Dmeasure = dsnormalize(Dmeasure)

% determine the number of tested features and the number of measures for
% each feature
[Nf,Nm] = size(Dmeasure);

% Compute the average measures
Dmean = mean(Dmeasure, 1); % row vector
Dmean = ones(Nf,1) * Dmean;

iDN = Dmeasure >= Dmean;

Dmeasure = 2 / Nm * (1- Dmeasure ./ Dmean);
Dmeasure(iDn) = 0;

end