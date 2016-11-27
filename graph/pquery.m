function dist = pquery(cost,x,y,grid)

[N,M] = size(cost);
if nargin<3
	grid = meshgrid(1:N,1:M);
end


%%%%%%%%%%%


%% load a synthetic cost surface image of size 512x512
n = 512;
[cost,W] = load_potential_map('mountain', n);
[N,M] = size(cost);
% display cost surface
figure, imagesc(cost), colormap gray, axis image

nsamples_list = round(linspace(20,300,10));
landmark = [];

for i=1:length(nsamples_list)
    nbr_landmarks = nsamples_list(i);
    
	% sample points using farthest geodesic sampling
	landmark = perform_farthest_point_sampling( W, landmark, nbr_landmarks-size(landmark,2), options );
	
	% display selected sample points
	Z = zeros(N,M);
	k = sub2ind(size(W),landmark(1,:)',landmark(2,:)');
	Z(k) = 1;
	figure, imagesc(Z), colormap gray, axis image
	
    % compute the associated triangulation
    [D,Z,Q] = perform_fast_marching(W, landmark);
	% display geodesic distance map
	figure, imagesc(D), colormap gray, axis image
	% display the associated Voronoi diagram
	figure, imagesc(Q), axis image
	
	%   DL(:,:,i) is the distance map to the ith landmark point.
	DL = zeros(N,M,nbr_landmarks);
	for i=1:nbr_landmarks
		DL(:,:,i) = perform_fast_marching(W, landmark(:,i));
	end

	[D1,Z] = compute_distance_landmark(start_points, DL, landmark, landmark_method);
	
	start_point=rand(1,2);

end;