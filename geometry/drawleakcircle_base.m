%% DRAWLEAKCIRCLE_BASE - Base function for DRAWLEAKCIRCLE. 
%
%% Syntax
%     [coord, A] = DRAWLEAKCIRCLE_BASE(R, method, gap, ngap, M);
%
%% See also 
% Related: 
% <DRAWLEAKCIRCLE.html |DRAWLEAKCIRCLE|>, 
% <DRAWCIRCLE_BASE.html |DRAWCIRCLE_BASE|>, 
% <DRAWGESTALT_BASE.html |DRAWGESTALT_BASE|>, 
% <BRESENHAMLINE_BASE.html |BRESENHAMLINE_BASE|>.
% Called:
% <../../statistics/html/UNIQUEUNSORT.html |UNIQUEUNSORT|>,
% <matlab:webpub(whichpath('POL2CART')) |POL2CART|>,
% <matlab:webpub(whichpath('CUMSUM')) |CUMSUM|>,
% <matlab:webpub(whichpath('MESHGRID')) |MESHGRID|>.

%% Function implementation
%--------------------------------------------------------------------------
function [coord, A] = drawleakcircle_base(R, method, gap, ngap, M)

center = [M M];

N = gausscircle(R)+1;

rho = repmat(R, [N 1]); % we don't take too much risk
theta =  transpose(linspace(0,2*pi,N));
[x,y] = pol2cart(theta, rho);
n = length(x);

x = round(x + repmat(center(1),[n 1])); % rounds toward 0
y = round(y + repmat(center(2),[n 1]));

pts = uniqueunsort([x y], 'rows', 'first');
x = pts(:,1);  y = pts(:,2);
n = length(x);

if gap>n
    error('drawleakcircle_base:inputerror', 'too small circle for input gap');
elseif ngap*gap>n
    warning('drawleakcircle_base:inputwarning', 'too small circle for input pair (gap,ngap)');
    ngap = 1;
end

range = round(n / ngap);

switch method
    
    case 'grid'
        mask = padarray(false(1,gap),[0 range-gap],true,'post');
        mask = repmat(mask,[1 ngap]);
        
        if length(mask)>=n
            mask = mask(1:n);
        else
            mask = padarray(mask,[0 n-length(mask)],true,'post');
        end
        
    case 'rand'
        mask = true(n,1);
        pos = 1 + round((range-gap-1)*rand(ngap,1)); % we avoid overlaps
        pos = repmat(pos(:),[1 gap]) + repmat(0:gap-1,[ngap 1]);
        % pos randomly distributed in [1,range-gap]
        offset = repmat([0 1:(ngap-1)]' * range, [1 gap]);
        pos = pos + offset;
        mask(pos(:)) = false;
end

x = x(mask);
y = y(mask);
coord = [x, y];

A = false(2*M+1,2*M+1);
A(sub2ind(size(A),x + (y-1)*(2*M+1))) = true;

end % end of drawleakcircle_base


%% Subfunction

%--------------------------------------------------------------------------
function N = gausscircle(r)
% count the number of lattice points N inside the boundary of a circle of 
% radius r 
%		N(r) = 1 + 4 \floor(r) + ...
%              4 \sum_{i=1}^\floor(r) \floor(\srqt(r^2 - i^2))
% http://mathworld.wolfram.com/GausssCircleProblem.html

nr = size(r,1);
fr = floor(r);
frmax = max(fr);

R = repmat(r,[1 frmax]);
FR = meshgrid(1:frmax,1:nr);
FR = cumsum(floor(sqrt(R.^2 - FR.^2)),2);

N = 1 + 4*fr + 4*FR((fr-1)*nr+(1:nr)');
end
