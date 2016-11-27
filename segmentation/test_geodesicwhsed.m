%test_aniwshed


%% Step 1: load image
rgb = imread('pears.png');
I = rgb2gray(rgb);
imshow(I)

%% Step 2: Use the Gradient Magnitude as the Segmentation Function

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
figure, imshow(gradmag,[]), title('Gradient magnitude (gradmag)')

L = watershed(gradmag);
Lrgb = label2rgb(L);
figure, imshow(Lrgb), 
title('Watershed transform of gradient magnitude (Lrgb)')

%% Step 3: Mark the Foreground Objects

se = strel('disk', 20);
Io = imopen(I, se);
figure, imshow(Io), title('Opening (Io)')

% compute the opening-by-reconstruction using imerode and imreconstruct
Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);
figure, imshow(Iobr), title('Opening-by-reconstruction (Iobr)')


% opening followed with a closing can remove the dark spots and stem 
% marks
% Compare a regular morphological closing with a closing-by-reconstruction
Ioc = imclose(Io, se);
figure, imshow(Ioc), title('Opening-closing (Ioc)')

% use imdilate followed by imreconstruct. Notice you must complement the
% image inputs and output of imreconstruct.
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
figure, imshow(Iobrcbr), title('Opening-closing by reconstruction (Iobrcbr)')

% calculate the regional maxima of Iobrcbr to obtain good foreground
% markers.
fgm = imregionalmax(Iobrcbr);
figure, imshow(fgm), 
title('Regional maxima of opening-closing by reconstruction (fgm)')

% superimpose the foreground marker image on the original image.
I2 = I; I2(fgm) = 255;
figure, imshow(I2), 
title('Regional maxima superimposed on original image (I2)')

% closing followed by an erosion
se2 = strel(ones(5,5));
fgm2 = imclose(fgm, se2);
fgm3 = imerode(fgm2, se2);

% use bwareaopen to removes all blobs that have fewer than a certain number 
% of pixels.
fgm4 = bwareaopen(fgm3, 20);
I3 = I;
I3(fgm4) = 255;
figure, imshow(I3)
title('Modified regional maxima superimposed on original image (fgm4)')


%% Step 4: Compute Background Markers

% thresholding operation.
bw = im2bw(Iobrcbr, graythresh(Iobrcbr));
figure, imshow(bw),
title('Thresholded opening-closing by reconstruction (bw)')

% we don't want the background markers to be too close to the edges of the
% objects we are trying to segment. We'll "thin" the background by computing 
% the "skeleton by influence zones", or SKIZ, of the foreground of bw. 
% compute the watershed transform of the distance transform of bw, and then 
% looking for the watershed ridge lines (DL == 0) of the result.
D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;
figure, imshow(bgm), title('Watershed ridge lines (bgm)')

%% Step 5: Compute the Watershed Transform of the Segmentation Function.

% use imimposemin  to modify an image so that it has regional minima only
% in certain desired locations. H
% modify the gradient magnitude image so that its only regional minima occur 
% at foreground and background marker pixels.
gradmag2 = imimposemin(gradmag, bgm | fgm4);

% compute the watershed-based segmentation.
L = watershed(gradmag2);


I4 = I;
I4(imdilate(L == 0, ones(3, 3)) | bgm | fgm4) = 255;
figure, imshow(I4)
title('Markers and object boundaries superimposed on original image (I4)')

Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
figure, imshow(Lrgb)
title('Colored watershed label matrix (Lrgb)')


%%%%%%%%%%%
rho = 0.5;
sigma = 1;
A = double(I);
T = gstsmooth(A,rho,sigma,'der','pey','sm','kro');

L = gstfeature(T(:,:,1,1), T(:,:,2,2), T(:,:,1,2),'norm','eign','sum');
figure, imshow(rescale(L,0,1)), title('tensor norm')
L = imimposemin(L, fgm4);

TT = 1 ./ (1+L);

% M = localextrema(L, 'min', 3);
% M = bwulterode(M);
% figure, imshow(M), title('Local minima of the tensor norm')

M = bwulterode(fgm4);

[i,j] = find(M);
% start_points must be of size 2 x nb_start_points.
start_points  = [i'; j'];

%start_points = [50 ;50];

[D,S,Q] = propagatefront(TT,start_points);


