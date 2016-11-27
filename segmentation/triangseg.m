function [vertex, faces, D, Q] = triangseg(I,varargin)

%% Parsing parameters

error(nargchk(1, 15, nargin, 'struct'));
error(nargoutchk(1, 5, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(I)
    error('gstsmooth:inputerror','a matrix is required in input'); 
end

% optional parameters
p = createParser('GSTFEATDETECT');   
% principal optional parameters
p.addOptional('sigma', 0.5, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addOptional('rho', 1, @(x)isscalar(x) && isfloat(x) && x>=0);
% additional optional parameters
p.addParamValue('der', 'pey', @(x)ischar(x) && ...
    any(strcmpi(x,{'mat','prew','cir','opt','sob',...
    'rad','dir','ana','pey','kov','kro'})));
p.addParamValue('sm', 'koe', @(x)ischar(x) && ...
    any(strcmpi(x,{'mat','pey','kro','iso','koe','ani'})));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%% Setting internal variables
[X,Y,C] = size(I);                                                     %#ok   
% possibly multichannel when C>1                            

A = expand2d(I,'cspl',3);

[GST,Dedge,Dcorner] = gstfeatdetect(A, p.sigma, p.rho,...
    'samp',1,'der',p.der,'sm',p.sm);
nGST = gstfeature(GST(:,:,1,1), GST(:,:,2,2), GST(:,:,1,2),...
    'norm','eign','dif');   

% Binary image of edges and corners
vertex = double(Dedge | Dcorner);
vertex = reduce2d(vertex,'cspl',3);

% Seeds points, include the corners into these points.
[i,j] = find(vertex);
vertex  = [i'; j'];
vertex(:,end+1:end+4) = [[1;1] [1;Y] [X;Y] [X;1]];

%faces = compute_delaunay(vertex);
%[vertex,D,Q] = perform_farthest_point_sampling( T, vertex, 10 );

alpha=2;
D = rescale(nGST,0,1-eps);
W = (1-D).^alpha;


options.null=0;

[D,S,Q] = perform_fast_marching(W, vertex, options);
Q = padarray(Q, [1 1], -1);
faces = compute_voronoi_triangulation(Q,vertex);
edges = compute_edges(faces);
ii = edges(1,:)>0 & edges(2,:)>0;
edges = edges(:,ii);

% plot sampling location
ms = 12; lw = 1.5; i = 0;

clf;
hold on;
imageplot(Q');
plot(edges(1,:), edges(2,:), 'r.', 'MarkerSize', ms);
hold off;
colormap jet(256);

clf;
hold on;
if not(isempty(edges))
    h = plot_edges(edges, vertex, 'b');
    set(h, 'LineWidth',lw);
end
plot(edges(1,:), edges(2,:), 'r.', 'MarkerSize', ms);
hold off;
axis tight; axis image; axis off;
colormap gray(256);

return
[D2,Z,Q2] = perform_fast_marching(W, landmark);

%[D,S,Q] = perform_fast_marching_mesh(vertex, faces, landmark, T);

% A way to compute the approximation is to compute coefficients vapprox that
% performs the best L2 approximation with linear spline.



end