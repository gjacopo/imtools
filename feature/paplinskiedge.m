function paplinskiedge(I)
% References:
%   [PB95]  A.P. Paplinski and J.F. Boyce: "Computational aspects of 
%       segmentation of a class of medical images using the concept of 
%       conjugate images", Tech. report 95-06, Monash University, 1995. 
%   [PK91]  M. Petrou and J. Kittler: "Optimal edge detectors for ramp
%       edges", IEEE Trans on Pattern Analysis and Machine Intelligence, 
%       13(5):483?491
%   [Spac86]  L.A. Spacek: "Edge detection and motion detection", Image &
%       Vision Computing, 4:43-56, 1986.

% Paplinski and Boyce used a two dimensional edge half operator h obtained 
% from a 1d filter introduced by Spacek [Spac86] and improved by Petrou [PK91]
% to perform the directional filtering in the four directions (0 90 180, 270)
% by rotating the filter mask. 
   
% In percentage values, the upper quarter of this filter is 
h = [ 0 8 11 8 2 0 0;...
0 16 26 26 13 1 0 ;...
0 26 48 52 33 8 0 ;...
0 43 80 80 53 17 0 ;...
0 68 100 94 61 21 0]; 

end