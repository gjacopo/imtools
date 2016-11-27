function [Q, B] =  ...
    ui_geodesicsegfgbg(I, method,  win, a, rho, sigma, der, int, samp, eign)
% UI_GEODESICSEGFGBG - Perform foreground/background segmentation using
% geodesic watershed.
%
%   [Q, B] =  ui_geodesicsegfgbg(I, method,  win, a, ...
%                                rho, sigma, der, int, samp, eign);
% See also 
% Related: 

% we accept a variable number of inputs
if nargin<10,     eign = 'l1';
    if nargin<9,     samp = 1;
        if nargin<8,     int = 'fast';
            if nargin<7,     der = 'fast';
                if nargin<6,     sigma = 1;
                    if nargin<5,     rho = 1;
                        if nargin<4,     a = [1 2];
                            if nargin<3,     win = 3;
                                if nargin<2,    method = 'ani'; end
                            end
                        end
                    end
                end
            end
        end
    end
end

C = size(I,3);

%% Define interactively the markers

disp('select 2 points in input image for foreground and background respectively');

figure, imagesc(rescale(I)), axis image, axis off
if C==1,  colormap gray;  end

markers = round(fliplr(ginput(2))'); 

close;


%% Perform geodesic watershed

disp('computing seeded geodesic watershed');
[D, Q] = geodesicwshed_base(I, markers, method, rho, sigma, a, der, int, samp, eign);


%% Define the background/foreground segments

sel = strel('square',3);
B = bwmorph(imdilate(Q,sel)-imerode(Q,sel)>0,'thin',Inf);

CC = repmat(B,[1 1 C]); 

IS = rescale(I) .* ~CC + CC;

figure, 
subplot(1,2,1), imagesc(Q), axis image, axis off, title('segments');
subplot(1,2,2), imagesc(IS), axis image, axis off, title('segmented regions');

end