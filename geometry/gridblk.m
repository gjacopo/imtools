%% GRIDBLK - Tiling of an image.
% 
%% Description
% Given an input image or the vector of its dimensions, compute the
% coordinates of the centers and the size of the square blocks dividing the
% image in a given number of regular tiles.
%
%% Syntax
%     [C,S] = GRIDBLK(I, nb);
%
%% Inputs
% *|I|* : image or vector of dimensions.
% 
% *|nb|* : upper bound for the total number of desired blocks dividing the
%     image regularly.
%
%% Outputs
% *|C|* : vector of coordinates of the centers of the tiling blocks.
% 
% *|S|* : vector of dimensions of the tiling blocks.
%
%% See also 
% Related: 
% <matlab:webpub(whichpath('BESTBLK')) |BESTBLK|>.
% Called:
% <matlab:webpub(whichpath('MESHGRID')) |MESHGRID|>.

%% Function implementation
function [C,S] = gridblk(I, nb)

XY = size(I(:,:,1));

if any(XY==ones(size(XY)))
    X = I(1);  Y = I(2);
else
    X = XY(1);  Y = XY(2);
end

S = sqrt(X* Y)/floor(sqrt(nb)); % nb will be an upper bound 
% S = sqrt(X* Y/nb); % nb will be a lower bound

%%
% compute the coordinates in x and y direction respectively
if true  % isempty(ver('images'))
    x = floor(X / S);
    if S*x<X,  x = x+1;  end
    if x==0,  x = floor(X/2);
    else      x = floor(S * (0:x-1) + S/2);
    end
    x(end) = min([X x(end)]);
    
    y = floor(Y / S);
    if S*y<Y,  y = y +1;  end
    if y==0,  y = floor(X/2);
    else      y = floor(S * (0:y-1) + S/2);
    end
    y(end) = min([Y y(end)]);
    S = [S S];

else
    sb = bestblk([X Y], S);                                            %#ok
    x = floor(S/2:sb(1):X);
    y = floor(S/2:sb(2):Y);
    S = [sb(1) sb(2)];
end

%%
% 'grid'
[x y] = meshgrid(x, y);

%%
% final vector of coordinates
C = [x(:) y(:)];

end % end of gridblk
