%% ELDERZUCKEREDGE_BASE - Base function for ELDERZUCKEREDGE.
%
%% Syntax
%     [edgemap, scmap, blur] = ELDERZUCKEREDGE_BASE(I, sigma);  
%
%% See also
% Related: 
% <ELDERZUCKEREDGE.html |ELDERZUCKEREDGE|>,
% <CANNYEDGE_BASE.html |CANNYEDGE_BASE|>,
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>,
% <CONGRUENCYEDGE_BASE.html |CONGRUENCYEDGE_BASE|>,
% <COMPASSEDGE_BASE.html |COMPASSEDGE_BASE|>,
% <ROTHWELLEDGE_BASE.html |ROTHWELLEDGE_BASE|>,
% <ANISOEDGE_BASE.html |ANISOEDGE_BASE|>,
% <KOETHEDGE_BASE.html |KOETHEDGE_BASE|>,
% <SDGDEDGE_BASE.html |SDGDEDGE_BASE|>,
% <PETROUEDGE_BASE.html |PETROUEDGE_BASE|>.
% Called: 
% <matlab:web(whichpath('EDGE')) |EDGE|>,
% <matlab:web(whichpath('CONV2')) |CONV2|>.

%% Function implementation
function [edgemap,varargout] = elderzuckeredge_base(I, sigma, reduce)  %#ok

%%
% dealing with multispectral images

C = size(I,3);
if C>1
    edgemap = zeros(size(I));
    if nargout>=2,  varargout{1} = zeros(size(I)); end
    if nargout==3,  varargout{2} = zeros(size(I)); end
    for c=1:C
        [edgemap(:,:,c),tmp1,tmp2] = elderzuckeredge_base(I(:,:,c), sigma);
        if nargout>=2,  varargout{1}(:,:,c) = tmp1; end
        if nargout>=3,  varargout{2}(:,:,c) = tmp2; end
    end
    return;
end

pad = 2;
I = padarray(I,[pad pad],'replicate','both');
[X,Y,C] = size(I);                                                     %#ok

%%
% computing edges                                                        

%alphaP = 1 - (1 - sigma).^(1./length(I))
alphaP = 2e-7;

mag = zeros(size(I));
ang = zeros(size(I));
map = zeros(size(I));

%%
% creating the different scales: these represent the different scale sizes
% supported here, expressed as the standard deviation of the filter.
sd1 = [16 8 4 2 1 0.5];

%%
% we iterate: for each standard deviation, we compute the 'gradient'- which
% in this case uses the steering functions defined by Elder and Zucker.
for h = 1:length(sd1)
    [mg ag map] = gradEZD2(I, sd1(h),sigma,alphaP, map);
    mag(mg ~= 0) = 0;
    mag = mag + mg;
    ang(ag ~= 0) = 0;
    ang = ang + ag;
end
%scaind = log2(map) + repmat(2, size(map));

sd2 = [4 2 1 0.5];

scmap = zeros([X,Y]);
lap_of_gau = zeros([X,Y]);

for h = 1:length(sd2)
    [l_o_g scmap] = lapEZD2(I,sd2(h),sigma,alphaP, ang, scmap);
    lap_of_gau(l_o_g ~= 0) = 0;
    lap_of_gau = lap_of_gau + l_o_g;
end
%scaind2 = log2(scmap) + repmat(2, size(scmap));

%%
% create the edgemap
edgemap = edge(lap_of_gau, 'canny');
edgemap = edgemap(pad+1:end-pad,pad+1:end-pad,:);

if nargout>=2
    varargout{1} = scmap(pad+1:end-pad,pad+1:end-pad,:);
end

%%
% create the blur map if required
if nargout==3
    
    % compute the distance between extrema
    [A,B] = find(edgemap ~= 0);
    [Alim, Blim] = size(I);
    window = 20;
    w = window;
    
    varargout{2} = zeros(X,Y);
    for i= 1:length(A),
        
        %define search window
        Alower = A(i) - window;
        Blower = B(i) - window;
        Aupper = A(i) + window;
        Bupper = B(i) + window;
        %clip to edge of image
        if Alower < 1,        Alower = 1;                              %#ok
        end
        if Blower < 1,        Blower = 1;                              %#ok
        end
        if Aupper>Alim,       Aupper = Alim;                           %#ok
        end
        if Bupper > Blim,     Bupper = Blim;                           %#ok
        end
        
        currentmax = lap_of_gau(A(i),B(i));
        currentmin = currentmax;
        nextmin = currentmin;
        nextmax = currentmax;
        
        currentpixelx = A(i);
        currentpixely = B(i);
        currentangle = 45*round(rad2deg(ang(currentpixelx, currentpixely))/45);
        if (currentangle == 0 )|| (currentangle == 360 )|| (currentangle == -360),
            pixelincx = 1;
            pixelincy = 0;
        end
        if (currentangle == 45 )|| (currentangle == -315),
            pixelincx = 1;
            pixelincy = 1;
        end
        if (currentangle == 90 )|| (currentangle == -270),
            pixelincx = 0;
            pixelincy = 1;
        end
        if (currentangle == 135 )|| (currentangle == -225),
            pixelincx = -1;
            pixelincy = 1;
        end
        if (currentangle == 180 )|| (currentangle == -180),
            pixelincx = -1;
            pixelincy = 0;
        end
        if (currentangle == 225 )|| (currentangle == -135),
            pixelincx = -1;
            pixelincy = -1;
        end
        if (currentangle == 270 )|| (currentangle == -90),
            pixelincx = 0;
            pixelincy = -1;
        end
        if (currentangle == 315 )|| (currentangle == -45),
            pixelincx = 1;
            pixelincy = -1;
        end
        counter = 0;
        while(nextmax >= currentmax),
            nextmaxx = currentpixelx + pixelincx;
            nextmaxy = currentpixely + pixelincy;
            if nextmaxx < 1 || (nextmaxx > size(I, 1)) || ...
                    (nextmaxy < 1) || (nextmaxy > size(I,2))
                break;
            end
            nextmax = lap_of_gau(nextmaxx, nextmaxy);
            counter = counter +1;
            if counter > w,
                break;
            end
            currentpixelx = nextmaxx;
            currentpixely = nextmaxy;
        end
        maxx = currentpixelx;
        maxy = currentpixely;
        
        currentpixelx = A(i);
        currentpixely = B(i);
        counter = 0;
        while(nextmin <= currentmin),
            nextminx = currentpixelx - pixelincx;
            nextminy = currentpixely - pixelincy;
            if nextminx < 1 || (nextminx > size(I, 1)) || ...
                    (nextminy < 1) || (nextminy > size(I,2))
                break;
            end
            nextmin = lap_of_gau(nextminx, nextminy);
            counter = counter +1;
            if counter > w,
                break
            end
            currentpixelx = nextminx;
            currentpixely = nextminy;
            
        end
        minx = currentpixelx;
        miny = currentpixely;
        
        d = sqrt((maxx-minx)^2 + (maxy-miny)^2);
        temp = sqrt((d/2)^2 - scmap(A(i),B(i))^2);
        if isreal(temp)
            varargout{2}(A(i),B(i)) = temp;
        else
            varargout{2}(A(i),B(i)) = 0;
        end
    end
    varargout{2} = varargout{2}(pad+1:end-pad,pad+1:end-pad,:);
   
end
end % end of elderzuckeredge_base


%% Subfunctions

%%
% |GRADEZD2| - Return two matrices the size of |im| that contain the magnitude
% and the angle of the intensity gradient of |im| using a gaussian directional
% derivative filter of standard deviation |sd1|.
% If given, the optional |marker| matrix will be set to |sd1| everywhere the
% gradient magnitude exceeds the critical value function of |sd1| and left 
% alone elsewhere.  This matrix will then be returned.  If no marker matrix
% is given as an argument, the function won't return the third value.
%--------------------------------------------------------------------------
function [mg, ag, marker, crit] = gradEZD2(im,scale,sigma,alphaP, marker)

if (nargin > 5)
  error('Wrong number of arguments to gradient.');
end
if (nargin > 4)
  if (size(marker) ~= size(im))
    error('Marker image not the same size as im.');
  end
else
  if (nargout > 2)
    error('No marker matrix can be returned unless you supply one.');
  end
end

% we want 2 sd on either side of the filter.
% this means that in all, our filter is 
%  4*scale by 4*scale
tail = ceil(2 * scale);

% construct the filter.
[X Y] = meshgrid(-tail:tail);
gx = g1x(X, Y, scale);
gy = g1x(Y, X, scale);
%convolve with the filter.
% this becomes very expensive with large sized filters.
gimx = conv2(im, gx, 'same');
gimy = conv2(im, gy, 'same');
ag = angle(gimx + 1i*gimy);
mg = steer1(ag, gimx, gimy);

abmag = abs(mg);
% compute the critical threshold for the gaussian of this scale.
% then threshold the magnitude by this value.
crit = c1(scale,sigma,alphaP);
if (nargin > 4)
  list = abmag >= crit;
  marker(list) = scale;
end

% output magnitude and angle of only those points 
%   with absolute magnitude of gaussian greater than
%   the critical threshold.
list = find(abmag < crit);
mg(list) = 0;
ag(list) = 0;
end % end of gradEZD2


%%
% |LAPEZD2| - Return a matrix representing the laplacian of |im| at the
% angles given by |gau_ang|. If given, the optional |marker| matrix will be 
% set to |sd2| everywhere the gradient magnitude exceeds the critical value
% function of |sd2| and left alone elsewhere.  This matrix will then be 
% returned. If no marker matrix is given as an argument, the function won't
% return the third value.
%--------------------------------------------------------------------------
function [lap_of_gau, marker] = lapEZD2(im,sd2, sigma,alphaP, gau_ang, marker)

if (nargin > 6)
  error('Wrong number of arguments to gradient.');
end
if (nargin > 5)
  if (size(marker) ~= size(im))
    error('Marker image not the same size as im.');
  end
else
  if (nargout > 1)
    error('No marker matrix can be returned unless you supply one.');
  end
end

% We want 2 sd on either side.
tail = ceil(2 * sd2);

[X Y] = meshgrid(-tail:tail);
gx = g2x(X, Y, sd2);
gy = g2x(Y, X, sd2);
gxy= g2xy(X, Y, sd2);
gimx = conv2(im, gx, 'same');
gimy = conv2(im, gy, 'same');
gimxy= conv2(im, gxy, 'same');
lap_of_gau = steer2(gau_ang, gimx, gimy, gimxy);

ablog = abs(lap_of_gau);
crit = c2(sd2,sigma,alphaP);
if (nargin > 5)
  list = ablog >= crit;
  marker(list) = sd2;
end

list = ablog < crit;
lap_of_gau(list) = 0;
end % end of lapEZD2


%%
% |STEER1|
%--------------------------------------------------------------------------
function imout = steer1(alpha, x_grad_im, y_grad_im)
imout = cos(alpha).*x_grad_im + sin(alpha).*y_grad_im;
end % end of steer1


%%
% |STEER2|
%--------------------------------------------------------------------------
function imout = steer2(alpha, xgrd_im, ygrd_im, xygrd_im)
ca = cos(alpha);
sa = sin(alpha);
imout =(ca.^2 .* xgrd_im)+(sa.^2 .* ygrd_im)-(2.*ca.*sa.*xygrd_im);
end % end of steer2


%%
% |G1X|
%--------------------------------------------------------------------------
function g = g1x(x,y,s1)
s1sq = s1.^2;
g = -(x./(2*pi*s1sq.^2)) .* exp(-(x.^2 + y.^2)./(2*s1sq)); 
end % end of g1x


%%
% |G2X|
%--------------------------------------------------------------------------
function g = g2x(x,y,s2)
sdnorm = 1 ./ (2*pi*s2.^4);
g = sdnorm .* (((x/s2).^2) - 1) .* exp(-(x.^2 + y.^2)./(2*s2.^2)); 
end % end of g2x


%%
% |G2XY|
%--------------------------------------------------------------------------
function g = g2xy(x,y,s2)
g = ((x.*y)./(2*pi*s2.^6)) .* exp(-(x.^2 + y.^2)./(2*s2.^2)); 
end % end of g2xy


%%
% |C1|
%--------------------------------------------------------------------------
function out = c1(sd1,sigma,alphaP)
s1 = sigma .* (1 ./ (2.*sqrt(2.*pi).*sd1.^2));
out = s1 .* sqrt(-2.*log(alphaP));
%out = 5 * out;
end % end of c1


%%
% |C2|
%--------------------------------------------------------------------------
function out = c2(sd2,sigma,alphaP)
s2 = sigma .* ((4.*sqrt(pi/3).*sd2.^3).^-1);
out = sqrt(2) .* s2 .* (erfinv(1-alphaP));
% best guess threshold
% out = (6 * sigma) ./ sd2.^3;
end % end of c2


%%
% |RAD2DEG|
%--------------------------------------------------------------------------
function degrees = rad2deg(radians)
degrees = radians * 180 / pi;
end % end of rad2deg
