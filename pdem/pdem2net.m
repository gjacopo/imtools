function [dem] = pdem2net(img,start_point,opt,imf,s1,s2,a)

[x,y] = size(img);
[X,Y] = meshgrid(1:x,1:y);
opt='i';
[dem,vdem] = pdem(img,start_point,opt,imf,s1,s2,a);

% Remove sinks with the function imfill
vdem = imfill(vdem,'holes');

% Upslope area (flow accumulation) algorithm for Digital Elevation Models
% calculate flow accumulation and direction
[flowacc, flowdir, slope, runs] = wflowacc(X,Y,dem,'type','multi');
% Derive Strahler Stream Order from Flow Direction Matrix
W = zeros(x,y);
W (start_point) = 1;
[S,nodes] = streamorder(flowdir,W);

subplot(1,2,1); 
pcolor(X,Y,+W); axis image; shading flat;
colorbar
title('Stream Network')
subplot(1,2,2);
pcolor(X,Y,S); axis image; shading flat;
colorbar
hold on


return

% test

% load (in this order)
addpath '/Users/236958/Work/MATLAB/'
addpath '/Users/236958/Work/MATLAB/peyre/toolbox_fast_marching/'
addpath '/Users/236958/Work/MATLAB/peyre/toolbox_fast_marching/toolbox/'
addpath '/Users/236958/Work/MATLAB/peyre/toolbox_diffc/'

imf=imread('/Users/236958/Work/Data/retina/retina-rank.tif'); imf=double(imf);
img=imread('/Users/236958/Work/Data/retina/retina-test.gif'); img=double(img);
[x,y] = size(img);
[X,Y] = meshgrid(1:x,1:y);
start_point=[225;120];
dem=pdem (img,start_point,'i',imf);
[flowacc, flowdir, slope, runs] = wflowacc(X,Y,dem,'type','single');
% Derive Strahler Stream Order from Flow Direction Matrix
W = zeros(x,y);
W (start_point(1,:),start_point(2,:)) = 1;
[S,nodes] = streamorder(flowdir,W);
