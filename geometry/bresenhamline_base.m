%% BRESENHAMLINE_BASE - Base function for BRESENHAMLINE. 
%
%% Syntax
%     [x y pts] = BRESENHAMLINE_BASE(x1, y1, x2, y2);
%
%% Acknowledgment
% Entirely derived from the implementation of A.Wetzler; see original source  
% code |BRESENHAM| included in this file, otherwise available at:
%   <http://www.mathworks.com/matlabcentral/fileexchange/28190>.
% Still, this original code did not support vectors of coordinates in input,
% contrary to the implementation |BRESENHAMLINE_BASE| proposed herein.
% See other implementations and source codes:
%   <http://www.mathworks.com/matlabcentral/fileexchange/12939>,
%   <http://www.mathworks.com/matlabcentral/fileexchange/13221>.
%
%% Credit
% <mailto:grazzja@lanl.gov J.Grazzini> (ISR-2/LANL)
%
%% See also 
% Related: 
% <BRESENHAMLINE.html |BRESENHAMLINE|>, 
% <SAMPLEDISK_BASE.html |SAMPLEDISK_BASE|>, 
% <SAMPLETRIANGLE_BASE.html |SAMPLETRIANGLE_BASE|>, 
% <DRAWLEAKCIRCLE_BASE.html |DRAWLEAKCIRCLE_BASE|>.
% Called:
% <matlab:webpub(whichpath('DIFF')) |DIFF|>,
% <matlab:webpub(whichpath('SUM')) |SUM|>,
% <matlab:webpub(whichpath('CUMSUM')) |CUMSUM|>,
% <matlab:webpub(whichpath('MESHGRID')) |MESHGRID|>.

%% Function implementation
%--------------------------------------------------------------------------
function [x, y, pts] = bresenhamline_base(x1, y1, x2, y2)

%%
% transform the inputs in column vector and round the coordinates
x1 = round(x1(:)); x2 = round(x2(:));
y1 = round(y1(:)); y2 = round(y2(:));

%%
% main computation starts here
s = size(x1(:),1);

dx = abs(x2-x1);
dy = abs(y2-y1);
steep = dy>dx; % vector column

% invert positions
t = dx(steep); dx(steep) = dy(steep); dy(steep) = t;  

% length of the longest line
mdx = max(dx)+1;
odx = ones(1,mdx);

%%
% transform steep
steep = steep(:, odx);   % repmat(steep, [1 mdx]);

I = dy==0;
DX = dx(:, odx);   % repmat(dx, [1 mdx]);
DY = dy(:, odx);   % repmat(dy, [1 mdx]);
A = floor(DX/2);
B = DX .* DY;
Z = cumsum(DY,2) - DY; 
Z = -Z;

Q = A+Z;
Q(Q < A-B) = 0;

Q = [zeros(s,1) diff(mod(Q,DX),1,2)>=0];
Q = cumsum(Q,2);
Q(I,:) = 0;

Z = meshgrid(0:mdx-1,1:s); % cumsum(ones(s, mdx),2)-1 

Y1 = y1(:, odx);   % repmat(y1, [1 mdx]);
X1 = x1(:, odx);   % repmat(x1, [1 mdx]);

Y2 = y2(:, odx);   % repmat(y2, [1 mdx]);
X2 = x2(:, odx);   % repmat(x2, [1 mdx]);

%%
% vectorized version of Bresnham algorith (see A.Wetzler's code)
A = X1<=X2;  pts = ~A;
x  = X1  + Z;
x(pts) = x(pts) - 2*Z(pts);
x(steep & A) = X1(steep & A) + Q(steep & A);
x(steep & pts) = X1(steep & pts) - Q(steep & pts); 
B = (A & x>X2) | (pts & x<X2);

A = Y1<=Y2;  pts = ~A;
y = Y1  + Z;
y(pts) = y(pts) - 2*Z(pts);
steep = ~steep;
y(steep & A) = Y1(steep & A) + Q(steep & A);
y(steep & pts) = Y1(steep & pts) - Q(steep & pts);
A = (A & y>Y2) | (pts & y<Y2);

%%
% trim
A = A | B;  
% additional output: but can be easily handled also without (using isnan)
pts = ~A; % pts = ~isnan(x);
x(A) = NaN;
y(A) = NaN;

end % end of bresenhamline_base


%% Original

%%
% |BRESENHAM| - Original implementation by A.Wetzler: single entries only
% (coordinates of two vertices to be linked by a line).
%--------------------------------------------------------------------------
function [x y] = bresenham(x1,y1,x2,y2)                                %#ok
%Matlab optmized version of Bresenham line algorithm. No loops.
%Format:
%               [x y]=bresenham(x1,y1,x2,y2)
%
%Input:
%               (x1,y1): Start position
%               (x2,y2): End position
%
%Output:
%               x y: the line coordinates from (x1,y1) to (x2,y2)
%
%Usage example:
%               [x y]=bresenham(1,1, 10,-5);
%               plot(x,y,'or');
x1=round(x1); x2=round(x2);
y1=round(y1); y2=round(y2);
dx = abs(x2-x1);
dy = abs(y2-y1);
steep = dy>dx;   % steep=abs(dy)>abs(dx)

if steep, t=dx; dx=dy; dy=t;  end

% the main algorithm goes here.
if dy==0 
    q=zeros(dx+1,1);
else
    q=[0;diff(mod((floor(dx/2):-dy:-dy*dx+floor(dx/2))',dx))>=0];
end

% and ends here.
if steep
    if y1<=y2,  y = (y1:y2)'; 
    else        y = (y1:-1:y2)';   end
    if x1<=x2,  x = x1+cumsum(q); 
    else        x = x1-cumsum(q);  end
else
    if x1<=x2,  x = (x1:x2)';  
    else        x = (x1:-1:x2)';   end
    if y1<=y2,  y = y1+cumsum(q);
    else        y = y1-cumsum(q);  end
end

end
