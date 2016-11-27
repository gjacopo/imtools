%% HIERTOPOSCALE_BASE - NOT IMPLEMENTED
%
%% Syntax
%   [MCC, S] = hiertoposcale_base(I, Hier, lambda)
%
%% Reference
%   [LAG09]  B. Luo,  J.-F. Aujol and Y. Gousseau: "Local scale measure from
%      the topographic map and application to remote sensing images", SIAM
%      Multiscale Modeling and Simulation, 8(1):1-29, 2009. 
%      http://www.cmi.univ-mrs.fr/~aujol/HDR/A16.pdf

%% Function implementation
function [MCC, S] = hiertoposcale_base(I, Hier, lambda)

[X,Y,C] = size(I);

nlevels = length(Hier);
% MC = cell(nlevels, 1);
MC = zeros(X, Y, nlevels);

% define the most contrasted component: this is done for the first level of
% the hiearchy, not pixel by pixel (in fact, all the pixels belonging to a
% same component in the first level will be assigned the same most contrasted)

% 'Area' : Scalar; the actual number of pixels in the region. (This value
%   might differ slightly from the value returned by bwarea, which weights
%   different patterns of pixels differently.)
% 'ConvexImage' : Binary image (logical) that specifies the convex hull, 
%   with all pixels within the hull filled in (i.e., set to on). (For pixels 
%   that the boundary of the hull passes through, regionprops uses the same
%   logic as roipoly to determine whether the pixel is inside or outside the 
%   hull.) The image is the size of the bounding box of the region. This 
%   property is supported only for 2-D input label matrices.
% 'ConvexArea' : Scalar that specifies the number of pixels in 'ConvexImage'. 
%   This property is supported only for 2-D input label matrices. 
% 'Perimeter' : p-element vector containing the distance around the boundary 
%   of each contiguous region in the image, where p is the number of regions. 
%   regionprops computes the perimeter by calculating the distance between
%   each adjoining pair of pixels around the border of the region. If the
%   image contains discontiguous regions, regionprops returns unexpected 
%   results. The following figure shows the pixels included in the perimeter
%   calculation for this object.
% 'PixelIdxList' : p-element vector containing the linear indices of the 
%   pixels in the region. 
% 'Solidity' : Scalar specifying the proportion of the pixels in the convex
%   hull that are also in the region. Computed as Area/ConvexArea. This
%   property is supported only for 2-D input label matrices

% first level
props = regionprops(Hier{1},{'Perimeter','Area','PixelIdxList'}); 

fi0 = props.Area;
per0 = props.Perimeter;

IL0 = I;
IL1 = zeros(size(I));

for l=2:nlevels
    
    props = regionprops(Hier{l},{'Perimeter','Area'});
    % areas of the componentes of the current level
    fi1 = props.Area;
    
    % spectral values of those components
    for i=1:C
        spec = accumarray(Hier{l}(:), reshape(I(:,:,c),[X*Y,1]), [], @mean);
        IL1(:,:,c) = reshape(spec, [X,Y]);
    end
    
    % test for the most contrasted: Eq.(3.2) of [LAG09]
    isContrasted = sign(fi1 - fi0 - lambda * per0) == -1;
    
    % contrast: Eq.(2.1)
    contrast = distcolor(IL1, IL0);
    % most contrasted: Eq.(3.3)
    MC(:,:,l) = contrast + isContrasted .* C{l-1};
    
    % update for the next level
    fi0 = fi1;
    IL0 = IL1;
    % perimeters (contours' lengths) used by the next level
    per0 = props.Perimeter;

end

[~, lMC] = max(MC, [], 3);
llevels = unique(lMC(:));


for l=ilevels(1):ilevels(end)
end
    



iCC = props.PixelIdxList;
% index = reshape(1:X*Y, [X,Y]);
% accumarray(MC, index, [], @(x){x});

end


function D = distcolor(A,B)
D = sqrt((A-B).^2);
% D = sum(abs(A-B),[],3);
end
