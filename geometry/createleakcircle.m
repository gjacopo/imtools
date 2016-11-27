function [x, y, A] = createleakcircle(R, method, gap, ngap, M, disp)
% CREATELEAKCIRCLE - Create a discrete disk (ie. with points on the regular 
% lattice) with approximated number of holes of approximated size. 
%
%    [x, y, A] = CREATELEAKCIRCLE(R, method, gap, ngap, M, disp);
%
% Inputs:
%   R : radius of the circle.
%   method : string defining the method used for creating holes in the set
%     of lattice points sampled on the disk; it is either:
%         - 'grid' for regular holes distribution,
%         - 'rand' for (pseudo) random holes distribution.
%   gap : size of the holes (in number of pixels/points).
%   ngap : desired number of holes.
%   M : half-size of the output image A displaying the circle; M should be
%     >(R+1); default: M=R+1.
%   disp : optional flag for drawing (discrete) lines; default: disp=false.
%
% Outputs:
%   x, y : (x,y) coordinates of the sampled points.
%   A : output logical image of size [2*M+1,2*M+1] with the 'leaking' circle
%     displayed (set to true for points whose coordinates are given by (x,y)).
% 
% credit: J.Grazzini (ISR-2/LANL)
%
% See also 
% Related: SAMPLEDISK, BRESENHAMLINE.

error(nargchk(4, 6, nargin, 'struct'));
error(nargoutchk(0, 3, nargout, 'struct'));

% set default
if nargin<6,  disp = false;  end

if nargin<5 || isempty(M)
    M = R+1;
elseif M<=R
    warning('createcircleleak:inputwarning', 'too small output image');
    M = R+1;
end

center = [M M];

N = gausscircle(R)+1;

rho = repmat(R, [N 1]); % we don't take too much risk
theta =  transpose(linspace(0,2*pi,N));
[x,y] = pol2cart(theta, rho);
n = length(x);

x = round(x + repmat(center(1),[n 1])); % rounds toward 0
y = round(y + repmat(center(2),[n 1]));

pts = unique([x y], 'rows');
x = pts(:,1);  y = pts(:,2);
n = length(x);

if gap>n
    error('createcircleleak:inputerror', 'too small circle for input gap');
elseif ngap*gap>n
    warning('createcircleleak:inputwarning', 'too small circle for input (gap,ngap)');
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

A = false(2*M+1,2*M+1);
A(sub2ind(size(A),x + (y-1)*(2*M+1))) = true;

if disp
    figure,   plot(x(:),y(:),'g.');
    axis([center(1)-1.1*R center(1)+1.1*R center(2)-1.1*R center(2)+1.1*R])
    figure, imagesc(A), colormap gray, axis image off
end

end


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
