%% SAMPLETRIANGLE - Sampling triangles over a grid.
% 
%% Description
% Sample lattice points (ie. with integer coordinates) inside triangles
% given by their respective faces and vertices coordinates.
%
%% Syntax
%    SAMPLETRIANGLE(tri, vertex);
%    [x,y] = SAMPLETRIANGLE(tri, vertex, method, samp, disp);
%
%% Inputs
% *|tri|* : matrix |(n,3)| storing the faces of |n| triangles; this is typically
%     the output of a triangulation.
% 
% *|vertex|* : matrix |(k,2)| storing the (X,Y) coordinates of the |k| vertices
%     indexed by the faces in |tri|.
% 
% *|method|* : (optional) string defining the method used for sampling the
%     lattice points inside the disks; it can be:
% 
% * |'full'| for all lattice points included in the triangle,
% * |'vista'| for original VISTA implementation,
% * |'grid'| for a variant of VISTA,
% * |'rand'| for a random selection of the lattice points;
% 
% note that in the cases |method='vista'|, |'grid'| or |'rand'|, the sample
%     points may lay outside of the triangles (small triangles), not in the
%     case |'full'|; default: |method='full'|.
% 
% *|samp|* : (optional) scalar used when |method='vista'|, |'grid'| or |'rand'|
%     (of no use when |method='full'|); it gives the number of sampled points 
%     when |samp>=1|, or the proportion of sampled points when |0<samp<1|;
%     in the latter case, the total number of points in a triangle is 'guessed'
%     using Pick's theorem; default: |samp=0.5|.
% 
% *|disp|* : (optional) flag for drawing (discrete) lines; default: |disp=false|.
% 
%% Outputs
% *|x, y : (X,Y) coordinates of the sampled points.
%
%% Reference
% See discussion at http://mathworld.wolfram.com/TrianglePointPicking.html
%
%% See also 
% Related: 
% <SAMPLEDISK.html |SAMPLEDISK|>, 
% <BRESENHAMLINE.html |BRESENHAMLINE|>.
% Called:
% <SAMPLETRIANGLE_BASE.html |SAMPLETRIANGLE_BASE|>.

%% Function implementation
%--------------------------------------------------------------------------
function varargout = sampletriangle(tri, vertex, method, samp, disp)

%%
% check arguments

error(nargchk(2, 5, nargin, 'struct'));
error(nargoutchk(0, 2, nargout, 'struct'));

% set default
if nargin<5  
    if nargout==0,  disp = true;
    else            disp = false;     end
    if nargin<4,  samp = [];
        if nargin<3,  method = 'full';  end
    end
end

if ~any(strcmpi(method,{'full','grid','vista','rand'}))
    error('sampletriangle:methoderror',['unknown method ' method]);

elseif isempty(samp) && any(method,{'grid','vista','rand'})
    % error('sampledisk:errorinput', ...
    %  'samp parameter must be entered with methods ''grid'', ''vista'' or ''rand''');
   samp = 0.5;       

elseif samp<=0 
    error('sampletriangle:errorinput', 'samp parameter must be >0');
end

%% 
% main computation

[x, y] = sampletriangle_base(tri, vertex, method, samp);

%%
% display

if disp
    plotsampletriangle(x, y, tri, vertex);
end

if nargin>=1,  varargout{1} = x;
    if nargin>=2,  varargout{2} = y;  end
end

end % end of sampletriangle


%% Subfunction

%--------------------------------------------------------------------------
function plotsampletriangle(x, y, tri, vertex)                        
ntri = size(tri,1);

for i=1:ntri
    figure % plot the triangles in red, the  sampled points as yellow dots
    plot(vertex(tri(i,[1 2 3 1]),1)',vertex(tri(i,[1 2 3 1]),2)','ro', ...
        vertex(tri(i,[1 2 3 1]),1)',vertex(tri(i,[1 2 3 1]),2)','r-', ... 
      x(i,:),y(i,:),'b.')
    % axis([]);
end
end