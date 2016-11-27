%% DRAWGESTALT_BASE - Base function for DRAWGESTALT.
% 
%% Syntax
%     DRAWGESTALT();
%     [contours, shapes] = DRAWGESTALT();
%     [contours, shapes] = DRAWGESTALT(X, Y);
%
%% See also 
% Related: 
% <DRAWGESTALT.html |DRAWGESTALT|>, 
% <DRAWCIRCLE_BASE.html |DRAWCIRCLE_BASE|>.
% <DRAWLEAKCIRCLE_BASE.html |DRAWLEAKCIRCLE_BASE|>.
% Called:
% <BRESENHAMLINE_BASE.html |BRESENHAMLINE_BASE|>.

%% Function implementation
function [C, S] = drawgestalt_base(X, Y, gap)

%%
% checking/setting parameters

error(nargchk(1, 3, nargin, 'struct'));
error(nargoutchk(0, 2, nargout, 'struct'));

% X=150; Y=200;
if nargin<3,  gap = round(X/4); 
    if nargin<2,  Y = [];  end
end

if isempty(Y),  Y = X;  end % default square size

C = struct('I',[], 'X',[], 'E',[], 'Y',[], 'T',[], 'L',[]);
S = struct('N',[], 'N2',[], 'F',[], 'Y',[], 'T',[], 'L',[]);

%% 
% draw incomplete contours

%%
% * junction type 1: endlink1
C.I = false(X,Y); 
C.I(1:round(X/2)-10,round(Y/2)) = true; 
C.I(X-20,10:Y-10) = true;

%%
% * intersection type L
C.L = false(X,Y);
C.L(20:round(X/2)+20,round(Y/2)-round(gap/2)) = true; 
C.L(round(X/2)+20,10:round(Y/2)-round(gap/2)) = true;
C.L(round(X/2)+20,round(Y/2)+round(gap/2):Y-10) = true;
C.L(X-10,10:Y-10) = true;

%%
% * intersection type Y
C.Y = false(X,Y);
[x,y] = bresenhamline([10 10 round(X/2)+20], ...
    [10 Y-10 round(Y/2)], ...
    [round(X/2)-20 round(X/2)-20 X-10], ...
    [round(Y/2)-20 round(Y/2)+20 round(Y/2)]);
C.Y(x+(y-1)*X) = true;

%%
% * intersection type T
C.T = false(X,Y);
C.T(40,10:round(Y/2)-round(gap/2)) = true;
C.T(40,round(Y/2)+round(gap/2):Y-10) = true;
C.T(40+round(gap/2):X-10,round(Y/2)) = true;

%%
% * multiple junction type 1
C.X = false(X,Y); 
C.X(1:round(X/2)-20,round(Y/2)) = true; 
C.X(round(X/2),10:round(Y/2)-20) = true;
C.X(round(X/2),round(Y/2)+20:Y-10) = true;
C.X(round(X/2)+20:end,round(Y/2)) = true; 

%%
% * incomplete contour
C.E1 = false(X,Y);
C.E1(20,10:Y-10) = true;
C.E1(X-20,10:Y-10) = true;
C.E1(round(X/2),10:round(Y/2)-round(gap/2)) = true;
C.E1(round(X/2),round(Y/2)+round(gap/2):Y-10) = true;

%%
% * composed incomplete contours
C.E2 = false(X,Y);
C.E2(10:round(X/2)-round(gap/2),round(Y/2)) = true; 
C.E2(round(X/2),10:round(Y/2)-round(gap/2)) = true;
[x,y] = bresenhamline(round(X/5), Y-10, ...
    round(X/2), round(Y/2)+round(gap/2));
C.E2(x(1,:)+(y(1,:)-1)*X) = true;
[x,y] = bresenhamline(round(X/2)+round(gap/2), round(Y/2),...
     X-10, round(Y/5));
C.E2(x+(y-1)*X) = true;

% closure
% [x, y, C] = createleakcircle(50, 'grid', 10, 10, 75);

%% 
% draw part shapes

%%
% * single-sided neck
S.N1 = false(X,Y);  
[x,y] = bresenhamline([10 X-10 round(X/2)-10 round(X/2)+10], ...
    [10 10 round(Y/2) round(Y/2)], ...
    [round(X/2)-10 round(X/2)+10 10 X-10], ...
    [round(Y/2) round(Y/2) Y-10 Y-10]);
S.N1(x(1:3,:)+(y(1:3,:)-1)*X) = true;
S.N1(round(X/2)+10, round(Y/2):Y-10) = true;


%%
% * double-sided neck
S.N2 = false(X,Y);
S.N2(x+(y-1)*X) = true;

%%
% * tapering shape
S.F = false(X,Y);
S.F(x(1:2,:)+(y(1:2,:)-1)*X) = true;
S.F(round(X/2)-10, round(Y/2):Y-10) = true;
S.F(round(X/2)+10, round(Y/2):Y-10) = true;

%%
% * intersection type Y
S.Y = false(X,Y);
S.Y(x([1,3],:)+(y([1,3],:)-1)*X) = true;
[x,y] = bresenhamline([round(X/2)-15 round(X/2)+25], ...
    [10 round(Y/2)+20], ...
    [round(X/2)+25 round(X/2)-15], ...
    [round(Y/2)-20 Y-10]);
S.Y(x+(y-1)*X) = true;
S.Y(round(X/2)+25:X-10,round(Y/2)-20) = true;
S.Y(round(X/2)+25:X-10,round(Y/2)+20) = true;

%%
% * intersection type T
S.T = false(X,Y);
% gap = 2*gap-20;
S.T(20,10:Y-10) = true;
S.T(round(X/2),10:round(Y/2)-round(gap/2)) = true;
S.T(round(X/2),round(Y/2)+round(gap/2):Y-10) = true;
S.T(round(X/2):X-10,round(Y/2)-round(gap/2)) = true;
S.T(round(X/2):X-10,round(Y/2)+round(gap/2)) = true;

%%
% * intersection type L
S.L = false(X,Y);
% gap = 2*gap-20;
S.L(20:X-10,20) = true;  
S.L(X-10,20:round(2*Y/3)) = true;
S.L(20:X-10-gap,20+gap) = true; 
S.L(X-10-gap,20+gap:round(2*Y/3)) = true;


end % end of drawgestalt