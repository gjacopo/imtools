%% BRESENHAMLINE - 2D Bresenham line algorithm.
% 
%% Description
% Apply the 2D Bresenham line algorithm [Bres65] to join 2 points of a grid.
% Given the list of lattice points (ie. points with integer coordinates) 
% linking two (or more) input vertices given by their (X,Y)-coordinates. 
%
% It is fully implemented in Matlab, optimized for vector manipulation (no 
% loops) and handles vectors of coordinates in input (contrary to other 
% original implementations).
%
%% Syntax
%     [x y pts] = BRESENHAMLINE(x1, y1, x2, y2);
%     [x y pts] = BRESENHAMLINE(x1, y1, x2, y2, disp);
%
%% Inputs
% *|(x1, y1)|* : start position(s) in the standard (X,Y)-orientated coordinate
%     system (X: horizontal OE and Y: vertical SN); can be vectors (row or
%     column) of vertices.
% 
% *|(x2, y2)|* : end position(s); must be of the same size as |(x1,y1)|.
% 
% *|disp|* : (optional) flag for drawing the estimated (discrete) lines joining
%     the |(x1,y1)| and |(x2,y2)| points; default: |disp=false|.
%
%% Outputs
% *|[x y]|* : the line coordinates from |(x1,y1)| to |(x2,y2)|; if the input
%      vertices are vectors, the list of coordinates for all lines is output;
%      each row of |x| (resp. |y|) provides the X(resp. Y)-coordinates of the
%      points on the line joining the vertices whose coordinates are given 
%      in the corresponding rows of |(x1,y1,x2,y2)|; the number of columns
%      of |x| (ibid |y|) is the maximum number of points found in a line; 
%      when a line has less points, it is complemented with |NaN| values.
% 
% *|pts|* : logical matrix of same size as |x| (and |y|) set to |true| for
%      the points really belonging to the line, |false| for the completing
%      |NaN| values: it is nothing else than |~isnan(x)|.
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
%% Example
%    x1 = [5 16 7 18 9];  y1 = [1 4 2 19 21]'; % whatever, row or column
%    x2 = [10 23 4 26 23];  y2 = [20 14 13 19 15]; % need to be of same size
%    bresenhamline(x1, y1, x2, y2, true);
%
%% References
% [Bres65]  J.E. Bresenham: "Algorithm for computer control of a digital
%      plotter", IBM Systems Journal, 4(1):25-30, 1965.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5388473&tag=1>
% 
% [wiki]  <http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm>
%
%% Credit
% <mailto:grazzja@lanl.gov J.Grazzini> (ISR-2/LANL)
%
%% See also 
% Related: 
% <SAMPLEDISK.html |SAMPLEDISK|>, 
% <SAMPLETRIANGLE.html |SAMPLETRIANGLE|>.
% Called:
% <BRESENHAMLINE_BASE.html |BRESENHAMLINE_BASE|>.

%% Function implementation
function [x, y, pts] = bresenhamline(x1, y1, x2, y2, disp)

%%
% settting/checking variables

% check errors in muber of input/output variables
error(nargchk(4, 5, nargin, 'struct'));
error(nargoutchk(0, 3, nargout, 'struct'));

if nargin<5,  disp = false;  end

if ~(isvector(x1) && isvector(y1) && isvector(x2) && isvector(y2))
    error('bresenhamline:inputerror','vectorial inputs required');
end

% transform the inputs in column vector and round the coordinates
x1 = round(x1(:)); x2 = round(x2(:));
y1 = round(y1(:)); y2 = round(y2(:));

if ~(isequal(size(x1),size(y1)) && isequal(size(y1),size(x2)) && ...
        isequal(size(x2),size(y2)))
    error('bresenhamline:inputerror','inputs with same dimension required');
elseif ~islogical(disp)
    error('bresenhamline:inputerror','logical flag required');
end

%% 
% main computation starts here

[x, y, pts] = bresenhamline_base(x1, y1, x2, y2);

%% 
% display

if disp
    nline = size(x,1);
    arbitrarymaxnline = 100;  % we are cheating here, but plot would be too slow...
    if nline>arbitrarymaxnline
        warning('bresenhamline:plotwarning', ...
            [num2str(arbitrarymaxnline) 'to be displayed only']);
        nline = arbitrarymaxnline;
    end
    nline = min([nline; size(x,1)]);
    x1 = x(1:nline,:)';  y1 = y(1:nline,:)';
    figure, hold on,  plot(x1, y1, '-b');
    [i,j] = find(diff([pts(1:nline,:) zeros(nline,1)],[],2));
    [~,l]=sort(i);  j = j(l);
    plot([x1(1,1:nline); x1(j'+(0:nline-1)*size(x1,1))], ...
        [y1(1,1:nline); y1(j'+(0:nline-1)*size(x1,1))],'or');
    axis([min(x1(:))-1 max(x1(:))+1 min(y1(:))-1 max(y1(:))+1]);
    hold off
end
 
end % end of bresenhamline

