function varargout = sflowinf(X,Y,dem,varargin)

% SFLOWINF - Upslope area algorithm for Digital Elevation Models using
% single flow direction dinf
%
%      [flowdir,slope] = sflowinf(X,Y,dem)
%      ... = sflowinf(X,Y,dem,src,'propertyname',propertyvalue,...)
% 
% Single flow direction algorithm (flow occurs only along the steepest 
% descent) that routes through flat terrain (not sinks). 
% Remove sinks with the function imfill that requires the image processing 
% toolbox.
%
% Input: 
%   X,Y       coordinate matrices created by meshgrid
%   dem       digital elevation model same size as X and Y
%
% Properties:
% propertyname               propertyvalues
%   'seed'           index(ices), ie. position(s) in the dem of the source(s) of the
%                    network, ie. starting cells; should be a 1xn vector for n sources
%                    (possibly use sub2ind to convert x-y coordinates to indices)
%                    'seed' is required when using the method for nondispersive flow
%                    direction estimation; otherwise the highest point in the dem is
%                    arbitrarly selected. 
%   'mode'            see wflowacc
%   'routeflats'      see wflowacc
%
% Output:
% flowdir   flowdirection 
% slope     slope (sparse matrix)
% runs      number of loops to route across flats (scalar) 
%
%
% Example from Paik:
%   dem = [ 93  99  96  95  94; ...
%           94  95  97  96  93 ; ...
%           91  98  100 97  95; ...
%           92  94  96  98  94; ...
%           89  91  90  92  93]


if nargin<3;
    error('wrong number of input arguments')
else
    siz = size(dem);
    if any(size(X) ~= size(Y)) || any((size(X) ~= siz));
        error('X, Y and dem must have same size')
    end
end

% general values
[NX,NY] = size(dem);
nrc = numel(dem);
nans = isnan(dem);

% ******************************************************************
% Main computation

% computes the flow direction and downslope
%   for all the pixels in
[M, S] = dem_flow(dem);
Mnans=isnan(M(:));

% M contains the pixel flow direction, in radians, for each pixel of dem.  
% Pixel flow direction is measured counter clockwise from the east-pointing 
% horizontal axis.  M is NaN for each pixel of DEM that has no downhill
% neighboUrs.
% S contains the downward slope (along the pixel flow direction) for each 
% pixel of dem.  Negative values indicate that the corresponding pixel of E 
% has no downhill neighbours.


directions ={'east', 'eastnortheast', 'northeast', 'northnortheast', ...
    'north',  'northnorthwest', 'northwest', 'westnorthwest', ...
    'west', 'westsouthwest', 'southwest', 'southsouthwest', ...
    'south','southsoutheast', 'southeast', 'eastsoutheast'};
% { 'north', 'northnortheast', 'northeast', 'eastnortheast', ...
%    'east', 'eastsoutheast',  'southeast', 'southsoutheast', ...
%    'south', 'southsouthwest', 'southwest', 'westsouthwest', ...
%    'west',  'westnorthwest',  'northwest', 'northnorthwest' };

% What are the linear index offsets corresponding to each of the neighbors?
offset.north          = -1;
offset.northnortheast = NX - 2;
offset.northeast      = NX - 1;
offset.eastnortheast  = 2*NX - 1;
offset.east           = NX;
offset.eastsoutheast  = 2*NX + 1;
offset.southeast      = NX + 1;
offset.southsoutheast = NX + 2;
offset.south          = 1;
offset.southsouthwest = -NX + 2;
offset.southwest      = -NX + 1;
offset.westsouthwest  = -2*NX + 1;
offset.west           = -NX;
offset.westnorthwest  = -2*NX - 1;
offset.northwest      = -NX - 1;
offset.northnorthwest = -NX - 2;

Ooffset=[];
for d=1:numel(directions)
    dir = directions{d};
    Ooffset = [Ooffset; offset.(dir)];
end

% zones range in [1,16]
Zones = ceil(8*M/pi + 16*(M<=0));    
Zones = Zones(:);

id = (1:nrc)';
id = id(~Mnans);
downslope = id + Ooffset(Zones(~Mnans));
ii = find(downslope<1 | downslope>nrc);
downslope(ii) = []; 
id(ii) = []; 

% sparse representation
M = sparse(id,downslope,1,nrc,nrc);

% ******************************************************************
% creating output

varargout{1} =  M;
if nargout > 1;
    varargout{2} = S;
end

% and this is it...
return;


