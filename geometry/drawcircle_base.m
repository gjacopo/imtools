%% DRAWCIRCLE_BASE - Base function for DRAWCIRCLE.
%
%% Syntax
%     [coord, A] = DRAWCIRCLE_BASE(center, radius, NOP);
%
%% See also 
% Related: 
% <DRAWCIRCLE.html |DRAWCIRCLE|>, 
% <DRAWGESTALT_BASE.html |DRAWGESTALT_BASE|>.
% <DRAWLEAKCIRCLE_BASE.html |DRAWLEAKCIRCLE_BASE|>.
% Called:
% <matlab:webpub(whichpath('POL2CART')) |POL2CART|>,
% <matlab:webpub(whichpath('LINSPACE')) |LINSPACE|>,
% <matlab:webpub(whichpath('PLOT')) |PLOT|>.

%% Function implementation
function [coord, A] = drawcircle_base(center, radius, NOP)

if nargin<3,  NOP = 50;  end

THETA = linspace(0, 2*pi, NOP);
RHO = ones(1,NOP) * radius;
[X,Y] = pol2cart(THETA, RHO);
X = X + center(1);  X = X(:);
Y = Y + center(2);  Y = Y(:);

coord = [X, Y];

X = ceil(abs(X)); Y = ceil(abs(Y));
M = max([X;Y])+1;
A = false(M,M);
A(sub2ind(size(A), X + (Y-1)*M)) = true;

end % end of drawcircle_base