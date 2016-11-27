%% CANNYEDGEMAP_BASE - Base function for CANNYEDGEMAP.
%
%% Syntax
%     map = CANNYEDGEMAP_BASE(gx, gy, der, mag, or, hyst, ratio);
%     [map, mag, or] = CANNYEDGEMAP_BASE(gx, gy, der, mag, or, hyst, ratio);
%
%% Acknowledgment
% This function is a copy/paste and crop of several proposed functions 
% implemented for thresholding, non maximum suppression and thinning of edge
% maps. Specifically: VISTA functions, Kovesi's library.
% 
%% See also
% Related:
% <CANNYEDGEMAP.html |CANNYEDGEMAP|>,
% <CANNYEDGE_BASE.html |CANNYEDGE_BASE|>,
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>,
% <ROTHWELLEDGE_BASE.html |ROTHWELLEDGE_BASE|>,
% <CONGRUENCYEDGE_BASE.html |CONGRUENCYEDGE_BASE|>,
% <COMPASSEDGE_BASE.html |COMPASSEDGE_BASE|>,
% <ANISOEDGE_BASE.html |ANISOEDGE_BASE|>,
% <ELDERZUCKEREDGE_BASE.html |ELDERZUCKEREDGE_BASE|>,
% <KOETHEDGE_BASE.html |KOETHEDGE_BASE|>,
% <SDGDEDGE_BASE.html |SDGDEDGE_BASE|>,
% <PETROUEDGE_BASE.html |PETROUEDGE_BASE|>.

%% Function implementation
function [map,varargout] = cannyedgemap_base(gx, gy, der, mag, or, hyst, ratio)

%%
% retrieve the gradient magnitude if not passed as an argument
C = size(gx,3);

if isempty(mag)
    mag = zeros(size(gx(:,:,1)));
    for c = 1:C
        mag = max(mag, sqrt((gx(:,:,c).*gx(:,:,c)) + (gy(:,:,c).*gy(:,:,c))));
    end
    magmax = max(mag(:));
    if magmax>0
        mag = mag / magmax;   % normalize
    end
end

%%
% retrieve the gradient orientation if not passed as an argument
if isempty(or)
    or = atan2(gy,gx);  
end

if isempty(ratio),   ratio = [1/3 0.08];  end;

if strcmpi(der,'vista'),    nbins = 128;
else                        nbins = 64;   end;


%% 
% deal with the case we have features for each channel, ie. channels are in
% fact processed separately
C = size(mag,3);
if C>1
    map = zeros(size(mag));
    if nargout>=2,  varargout{1} = zeros(size(mag));
        if nargout==3,  varargout{2} = zeros(size(mag));  end
    end
    for c = 1:C
        [map(:,:,c), m, o] = cannyedgemap_base(gx(:,:,c), gy(:,:,c), der, ...
            mag(:,:,c), or(:,:,c), hyst, ratio);
    end
    if nargout>=2,  varargout{1}(:,:,c) = m;
        if nargout==3,  varargout{2}(:,:,c) = o;  end
    end
    return
end

%%
% compute hystheresis thresholds
if isempty(hyst)
    hyst = selecthyst(mag, ratio, nbins);
elseif hyst(1) > hyst(2)    % swap values
    tmp = hyst(1); hyst(1) = hyst(2); hyst(2) = tmp;
end

%%
% estimations

switch der
    
    case 'kovesi'
        radius = 1;
        % note that the orientation in Kovesi's function is defined with
        % respect to the standard Matlab axis
        or = (or+pi/2);  or(or>pi) = or(or>pi) - 2*pi;
        or = -or;        
        
        or = or.*~(or<0) + (or+pi).*(or<0); % map angles to 0-pi
        or = or * 180/pi;   % do convert to degrees
        map = nonmaxsup_kovesi(mag, or, radius);
        % hyst = max(mag(:)) * hyst; % mag already normalized
        map = hysthresh_kovesi(map, hyst);
        
    case {'matlab','vista'}
        map = hystnonmaxsup_matlab(gx, gy, mag, hyst );
end

%%
% outputs

if nargout>=2,  varargout{1} = mag;  
    if nargout==3,  varargout{2} = or;  end
end

end  % end of cannyedgemap_base


%% Subfunctions

%%
% |SELECTHYST| - Hysteresis thresholds' selection.
%--------------------------------------------------------------------------
function hyst = selecthyst(mag, ratio, nbins)
[m,n] = size(mag(:,:,1));
[counts,x] = imhist(mag, nbins);
pixratio = ratio(1);
if pixratio<0
    meanx = sum(counts.*x)/sum(counts);
    stdx = sqrt(sum(counts.*((x-meanx).^2))/sum(counts));
    pixratio = 1-(meanx+2*stdx)/max(x);
end
thresratio = ratio(2);
highThresh = find(cumsum(counts) > (1-pixratio)*m*n,1,'first') / nbins;
lowThresh = thresratio*highThresh;
hyst = [lowThresh highThresh];
end % end of selecthyst


%%
% |HYSTNONMAXSUP_MATLAB| - Non-maximum supression.
%--------------------------------------------------------------------------
function map = hystnonmaxsup_matlab(gx, gy, mag, hyst)
[m,n] = size(gx(:,:,1));                                               %#ok          
map = zeros(size(mag)); % map = repmat(false, m, n);

% we accrue indices which specify ON pixels in strong map; the array e will 
% become the weak edge map
idxStrong = [];
for dir = 1:4
    idxLocalMax = findmaxima_matlab(dir,gx,gy,mag);
    idxWeak = idxLocalMax(mag(idxLocalMax) > hyst(1));
    map(idxWeak)=1;
    idxStrong = [idxStrong; idxWeak(mag(idxWeak) > hyst(2))];          %#ok
end

rstrong = rem(idxStrong-1, m)+1;
cstrong = floor((idxStrong-1)/m)+1;
% apply connectivity rules
map = bwselect(map, cstrong, rstrong, 8);
% thin double (or triple) pixel wide contours
map = bwmorph(map, 'thin', 1);  
end % end of hystnonmaxsup_edge


%%
% |FINDMAXIMA_MATLAB| - This function helps with the non-maximum supression in 
% the Canny edge detector. It is a copy paste of the original FINDMAXIMA
% function implemented in Matlab.
% 
% Inputs:
%   |direction| : the index of which direction the gradient is pointing, read
%      from the diagram below. direction is 1, 2, 3, or 4.
%
%   |ix| : input image filtered by derivative of Gaussian along X
%
%   |iy| : input image filtered by derivative of Gaussian along Y
%
%   |mag| : the gradient magnitude image
%
% Output:
%   |idxLocalMax| : index of the local maxima in the gradient magnituce image
%     map
%--------------------------------------------------------------------------
function idxLocalMax = findmaxima_matlab(direction,ix,iy,mag)

[m,n] = size(mag(:,:,1));

%%
% find the indices of all points whose gradient (specified by the
% vector (ix,iy)) is going in the direction we're looking at.
switch direction
    case 1
        idx = find((iy<=0 & ix>-iy)  | (iy>=0 & ix<-iy));
    case 2
        idx = find((ix>0 & -iy>=ix)  | (ix<0 & -iy<=ix));
    case 3
        idx = find((ix<=0 & ix>iy) | (ix>=0 & ix<iy));
    case 4
        idx = find((iy<0 & ix<=iy) | (iy>0 & ix>=iy));
end

%%
% exclude the exterior pixels
if ~isempty(idx)
    v = mod(idx,m);
    idx(v==1 | v==0 | idx<=m | (idx>(n-1)*m)) = [];
end

ixv = ix(idx);
iyv = iy(idx);
gradmag = mag(idx);

%%
% do the linear interpolations for the interior pixels
switch direction
    case 1
        d = abs(iyv./ixv);
        gradmag1 = mag(idx+m).*(1-d) + mag(idx+m-1).*d;
        gradmag2 = mag(idx-m).*(1-d) + mag(idx-m+1).*d;
    case 2
        d = abs(ixv./iyv);
        gradmag1 = mag(idx-1).*(1-d) + mag(idx+m-1).*d;
        gradmag2 = mag(idx+1).*(1-d) + mag(idx-m+1).*d;
    case 3
        d = abs(ixv./iyv);
        gradmag1 = mag(idx-1).*(1-d) + mag(idx-m-1).*d;
        gradmag2 = mag(idx+1).*(1-d) + mag(idx+m+1).*d;
    case 4
        d = abs(iyv./ixv);
        gradmag1 = mag(idx-m).*(1-d) + mag(idx-m-1).*d;
        gradmag2 = mag(idx+m).*(1-d) + mag(idx+m+1).*d;
end
idxLocalMax = idx(gradmag>=gradmag1 & gradmag>=gradmag2);
end % end of findmaxima_edge


%%
% |HYSTHRESH_KOVESI| - Perform hysteresis thresholding of an image.
%
% Inputs:
%   |im| : image to be thresholded (assumed to be non-negative)
%
%   |hyst=[T1,T2]| : upper and lower threshold values |(T1 > T2)|
%
% Output:
%   |bw| : the thresholded image (containing values 0 or 1)
%--------------------------------------------------------------------------
function bw = hysthresh_kovesi(im, hyst)

T1 = hyst(1); T2=hyst(2);
if T1 < T2    % T1 and T2 reversed - swap values
    tmp = T1;
    T1 = T2;
    T2 = tmp;
end

aboveT2 = im > T2;   % edge points above lower threshold.
[aboveT1r, aboveT1c] = find(im > T1);  % Row and colum coords of points
% above upper threshold.

% Obtain all connected regions in aboveT2 that include a point that has a
% value above T1
bw = bwselect(aboveT2, aboveT1c, aboveT1r, 8);
end % end of hysthresh_kovesi


%%
% |NONMAXSUP_KOVESI| - Function for performing non-maxima suppression on an 
% image using an orientation image.  It is assumed that the orientation 
% image gives feature normal orientation angles in degrees (0-180).
%
% Inputs:
%   |inimage| : image to be non-maxima suppressed.
%
%   |orient| : image containing feature normal orientation angles in degrees
%      (0-180), angles positive anti-clockwise.
%
%   |radius| : distance in pixel units to be looked at on each side of each
%      pixel when determining whether it is a local maxima or not; this
%      value cannot be less than 1; suggested value about 1.2 - 1.5.
%
% Output
%   |im| : non maximally suppressed image.
%--------------------------------------------------------------------------
function im = nonmaxsup_kovesi(inimage, orient, radius)

[rows,cols] = size(inimage);
im = zeros(rows,cols); % preallocate memory for output image

iradius = ceil(radius);

%%
% precalculate x and y offsets relative to centre pixel for each orientation angle 

angle = (0:180).*pi/180; % array of angles in 1 degree increments (but in radians).
xoff = radius*cos(angle); % x and y offset of points at specified radius and angle
yoff = radius*sin(angle); % from each reference position.

hfrac = xoff - floor(xoff); % fractional offset of xoff relative to integer location
vfrac = yoff - floor(yoff); % fractional offset of yoff relative to integer location

orient = fix(orient)+1; % orientations start at 0 degrees but arrays start with index 1.

%%
% now run through the image interpolating grey values on each side
% of the centre pixel to be used for the non-maximal suppression.

for row = (iradius+1):(rows - iradius)
    for col = (iradius+1):(cols - iradius)
        
        or = orient(row,col);   % Index into precomputed arrays
        x = col + xoff(or);     % x, y location on one side of the point in question
        y = row - yoff(or);
        
        fx = floor(x);          % Get integer pixel locations that surround location x,y
        cx = ceil(x);
        fy = floor(y);
        cy = ceil(y);
        tl = inimage(fy,fx);    % Value at top left integer pixel location.
        tr = inimage(fy,cx);    % top right
        bl = inimage(cy,fx);    % bottom left
        br = inimage(cy,cx);    % bottom right
        
        upperavg = tl + hfrac(or) * (tr - tl);  % Now use bilinear interpolation to
        loweravg = bl + hfrac(or) * (br - bl);  % estimate value at x,y
        v1 = upperavg + vfrac(or) * (loweravg - upperavg);
        
        if inimage(row, col) > v1 % We need to check the value on the other side...
            
            x = col - xoff(or);     % x, y location on the `other side' of the point in question
            y = row + yoff(or);
            
            fx = floor(x);
            cx = ceil(x);
            fy = floor(y);
            cy = ceil(y);
            tl = inimage(fy,fx);    % Value at top left integer pixel location.
            tr = inimage(fy,cx);    % top right
            bl = inimage(cy,fx);    % bottom left
            br = inimage(cy,cx);    % bottom right
            upperavg = tl + hfrac(or) * (tr - tl);
            loweravg = bl + hfrac(or) * (br - bl);
            v2 = upperavg + vfrac(or) * (loweravg - upperavg);
            
            if inimage(row,col) > v2            % This is a local maximum.
                im(row, col) = inimage(row, col); % Record value in the output
                % image.
            end
            
        end
    end
end

%%
% finally thin the 'nonmaximally suppressed' image by pointwise
% multiplying itself with a morphological skeletonization of itself.
skel = bwmorph(im,'skel',Inf);
im = im.*skel;

end  % end of nonmaxsup_kovesi
