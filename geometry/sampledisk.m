%% SAMPLEDISK - Sampling disks over a grid.
% 
%% Description
% Sample lattice points (ie. with integer coordinates) inside disks given by
% their respective centers and radii.
%
%% Syntax
%    SAMPLEDISK(center, radius);
%    [x, y, w] = SAMPLEDISK(center, radius, method, samp, disp);
%
%% Inputs
% *|center|* : matrix |(k,2)| of the (X,Y) integer coordinates of the centers
%     of |k| disks.
% 
% *|radius|* : matrix |(k,1)| of the radii of the |k| disks.
% 
% *|method|* : (optional) string defining the method used for sampling the
%     lattice points inside the disks; it can be:
% 
% * |'full'| for all lattice points included in the disks,
% * |'grid'| for all lattices points regularly sampled on a 2D grid
%           with size equal to samp (see below) and included in the disks,
% * |'radi'| for points sampled on internal circles inside the disks
%          (by default, at half the size of the disk, otherwise depending 
%           on samp) plus the center of the circle itself,
% * |'rand'| for randomly selected points inside the disks;
% 
%     default: |method='full'|.
% 
% *|samp|* : (optional) scalar or vector of size |(k,1)| giving (of no use when
%     |method='full'| or |'scale'|):
%
% * for |method='grid'|, the sampling rate of the grid when |samp>1| or
%          the proportion of desired sampled points on this grid when 
%          |0<=samp<1|, 
% * for |method='rand'|, the total number of sampled points when |samp>=1|
%          and the proportion of sampled points when |0<=samp<1|,
% * for |method='radi'|, _(i)_ |when samp>1|, the total number of points
%          regularly sampled on the border of the disk, _(ii)_ when |samp<-1|, 
%          the total number of sampled points around an internal circle at
%          a distance radius/2 from the disk's center, _(iii)_ when |0<samp<=1|,
%          the ratio of the radius of an internal circle on which as many 
%          lattice points as possible are sampled; in particular, points are
%          sampled around the disk border when |samp=1|;
% 
% note that in the latter case (|'radi'|), the total number of points in the
%     disk is estimated using Gauss formula; default: |samp=0.5| for all
%     chosen methods. 
% 
% *|disp|* : (optional) flag for drawing (discrete) lines; default: |disp=false|.
% 
%% Outputs
% *|x, y|* : (X,Y) coordinates of the sampled points.
%
%% Reference
% See discussions at     http://mathworld.wolfram.com/DiskPointPicking.html
%                     http://mathworld.wolfram.com/GausssCircleProblem.html   
%
%% See also 
% Related: 
% <SAMPLETRIANGLE.html |SAMPLETRIANGLE|>, 
% <BRESENHAMLINE.html |BRESENHAMLINE|>.
% Called:
% <SAMPLEDISK_BASE.html |SAMPLEDISK_BASE|>.

%% Function implementation
%--------------------------------------------------------------------------
function varargout = sampledisk(center, radius, method, samp, disp)

%%
% settting/checking variables

% check errors in muber of input/output variables
error(nargchk(3, 5, nargin, 'struct'));
error(nargoutchk(0, 3, nargout, 'struct'));

% set default
if nargin<5  
    if nargout==0,  disp = true;
    else            disp = false;     end
    if nargin<4,  samp = [];
        if nargin<3,  method = 'full';  end
    end
end

if ~any(strcmpi(method,{'full','grid','scale','radi','rand'}))
    error('sampledisk:methoderror',['unknown method ' method]);
end

if isempty(samp) && any(strcmpi(method,{'grid','radi','rand'}))
    % error('sampledisk:errorinput', ...
    %  'samp parameter must be entered with methods ''grid'', ''radi'' or ''rand''');
   samp = 0.5;       
elseif ~isempty(samp) && any(strcmpi(method,{'full','scale'}))
    warning('sampledisk:warninginput', ...
        ['samp parameter ignored with method ' method]);
end

%% 
% main computation

[x, y, w] = sampledisk_base(center, radius, method, samp);

if disp
    plotsampledisk(x, y, center, radius);
end

if nargin>=1,  varargout{1} = x;
    if nargin>=2,  varargout{2} = y;
        if nargin==3,  varargout{3} = w;  end
    end
end

end % end of sampledisk


%% Subfunction

%--------------------------------------------------------------------------
function plotsampledisk(x, y, center, radius)                          
ndisks = size(center,1);

for i=1:ndisks
    figure
    plot(x(i,:),y(i,:),'g.')
    axis([center(i,1)-1.1*radius(i) center(i,1)+1.1*radius(i) ...
        center(i,2)-1.1*radius(i) center(i,2)+1.1*radius(i)])
    hold on
    t = linspace(0,2*pi,1000);  r = ones(1,1000)*radius(i);
    [xx,yy] = pol2cart(t,r);
    xx=xx+center(i,1);   yy=yy+center(i,2);
    plot(xx,yy,'--');
    hold off
    axis square
end
end
