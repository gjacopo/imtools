function Iint = normconv(Isamp, Csamp, g)
% Isamp: sampled image 
% g: smoothing ?lter 
% Csamp: sampling mask so that Isamp = I * Csamp

% step 1: convolve Isamp with g 
C = conv2(Isamp,g,'same');

% step (ii): convolve Csamp with g 
NC = conv2(Csamp,g,'same');

% step (iii): Divide the two results point by point: 
Iint = C ./ NC;

end