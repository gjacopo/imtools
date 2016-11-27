%% SAMPLETRIANGLE_BASE - Base function for SAMPLETRIANGLE.
%
%% Syntax
%    [x,y] = SAMPLETRIANGLE_BASE(tri, vertex, method, samp);
%
%% Reference
% See discussion at http://mathworld.wolfram.com/TrianglePointPicking.html
%
%% See also 
% Related: 
% <SAMPLETRIANGLE.html |SAMPLETRIANGLE|>, 
% <SAMPLEDISK_BASE.html |SAMPLEDISK_BASE|>, 
% <BRESENHAMLINE_BASE.html |BRESENHAMLINE_BASE|>.
% Called:
% <matlab:webpub(whichpath('POL2CART')) |POL2CART|>,
% <matlab:webpub(whichpath('MESHGRID')) |MESHGRID|>.

%% Function implementation
%--------------------------------------------------------------------------
function [x, y] = sampletriangle_base(tri, vertex, method, samp)

global include_vertex;
include_vertex = false;  % if true, then vertex will be included in the samples

if ~isempty(samp) && samp>=0 && samp<1,
    samp = floor(samp .* boundtripoints(tri, vertex));
end

switch method
    
    case 'full'
        [x, y] = fullsampletriangle(tri, vertex);

    case 'grid'
        samp = max(ceil((samp+2-include_vertex)/3),2);
        % we will found the number 3*n closest to samp
        [x, y] = gridsampletriangle(tri, vertex, samp);
        
    case 'vista'
        samp = samp + (1-include_vertex);
        [x, y] = vistasampletriangle(tri, vertex, samp);
    
    case 'rand'
        [x, y] = randsampletriangle(tri, vertex, samp);
                
    otherwise
        error('sampletriangle_base:methoderror',['unknown method ' method]);
end

end % end of sampletriangle_base


%% Subfunctions

%%
% |GRIDSAMPLETRIANGLE| - Regular sampling of a triangle over the underlying
% lattice.
%--------------------------------------------------------------------------
function [x, y] = gridsampletriangle(tri, vertex, samp)
global include_vertex

ntri = size(tri,1);
R = ((1-include_vertex)/samp:(1/samp):(1-(1/samp)));
nsamp = length(R); % depending on include_vertex, it is samp or (samp-1)
R = [R; R.*(1-R); 1-(2*R-R.^2)]';

X = reshape(vertex(tri(:,[ 1 2 3 2 3 1 3 1 2])',1),[3 3*ntri]);
Y = reshape(vertex(tri(:,[ 1 2 3 2 3 1 3 1 2])',2),[3 3*ntri]);

x = fix(transpose(reshape(R*X, [3*nsamp ntri ])));
y = fix(transpose(reshape(R*Y, [3*nsamp ntri ])));
% final number of points: 3*(samp-1) or 3*samp, depending on include_vertex

end


%%
% |VISTASAMPLETRIANGLE| - Sparse triangle sampling.
%--------------------------------------------------------------------------
function [x, y] = vistasampletriangle(tri, vertex, samp)
global include_vertex

ntri = size(tri,1);
R = ((1-include_vertex)/samp:(1/samp):(1-(1/samp)));
R = [R; R.*(1-R); 1-(2*R-R.^2)];

x = fix(reshape(vertex(tri,1),[ntri 3]) * R);
y = fix(reshape(vertex(tri,2),[ntri 3]) * R);

end


%%
% |VISTASAMPLETRIANGLE_NEW| - Modified sparse triangle sampling.
%--------------------------------------------------------------------------
function [x, y] = vistasampletriangle_new(tri, vertex, samp)           %#ok
global include_vertex

ntri = size(tri,1);

originx = vertex(tri(:,1),:);
originy = [zeros(ntri,1) originx(:,2) originx(:,2)];
originx = [zeros(ntri,1) originx(:,1) originx(:,1)];

R = ((1-include_vertex)/samp:(1/samp):(1-(1/samp)));
T =  fliplr(R);  % R + T <1
R = [ ones(1,length(R)); R; T ];

x = fix((reshape(vertex(tri,1),[ntri 3])-originx) * R);
y = fix((reshape(vertex(tri,2),[ntri 3])-originy) * R);

end


%%
% |FULLSAMPLETRIANGLE| - Full triangle sampling.
%--------------------------------------------------------------------------
function [x, y] = fullsampletriangle(tri, vertex)           
ntri = size(tri,1);

b = boundtripoints(tri, vertex); 
x = NaN(ntri, b);
y = NaN(ntri, b);

nsamp = 0;
for i=1:ntri
    xv = vertex(tri(i,:),1);
    yv = vertex(tri(i,:),2);
    xp = min(xv):1:max(xv);
    yp = min(yv):1:max(yv);
    [xg, yg] = meshgrid(xp, yp);
    ip = inpolygon(xg(:), yg(:), xv, yv);
    lip = sum(ip);
    x(i,1:lip) = xg(ip);
    y(i,1:lip) = yg(ip);
    nsamp = max([nsamp lip]);
end

x = x(:,1:nsamp);
y = y(:,1:nsamp);

end


%%
% |RANDSAMPLETRIANGLE| - Random triangle sampling.
%--------------------------------------------------------------------------
function [x, y] = randsampletriangle(tri, vertex, samp)
ntri = size(tri,1);

t = sqrt(rand(1,samp));
s = rand(1,samp);

R = [(1-s); s.*(1-t); s.*t];

x = fix(reshape(vertex(tri,1),[ntri 3]) * R);
y = fix(reshape(vertex(tri,2),[ntri 3]) * R);
end


%%
% |BOUNDTRIPOINTS| - Find a maximum possible size of a triangle sampling
% using Pick's theorem
%
% When the vertices of a triangle (true for any polygon in fact) are at 
% integer coordinates (lattice points) on a grid, then:
%
%      area + 1 = #{points inside triangle} + #{points on edge} / 2
%
% therefore, we are ensured that: 
%
%      #{points inside triangle} + #{points on edge} <= 2 * (area+1) 
%
% See also http://www.btinternet.com/~se16/hgb/triangle.htm
%--------------------------------------------------------------------------
function b = boundtripoints(tri, vertex)                       
b = max(abs( (vertex(tri(:,2),1) - vertex(tri(:,1),1)) .* ...
    (vertex(tri(:,3),2) - vertex(tri(:,1),2)) ...
    - (vertex(tri(:,2),2) - vertex(tri(:,1),2)) .* ...
    (vertex(tri(:,3),1) - vertex(tri(:,1),1)))) + 2;
end

