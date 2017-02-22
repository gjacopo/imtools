function [M,varargout] = sflow8(dem,varargin)
% SFLOW8 - Upslope area algorithm for Digital Elevation Models using
% single flow direction d8 and gd8
%
%     dir = sflow8(dem);
%     [dir,slope,upslope] = sflow8(dem, flow, ...
%                                  propertyname',propertyvalue,...);
% 
% Single flow direction algorithm (flow occurs only along the steepest
% descent) that routes through flat terrain (not sinks). 
% Remove sinks with the function imfill that requires the image processing 
% toolbox.
%
% Input: 
%   dem : digital elevation model same size as X and Y
%   flow : optional string setting the nature of the estimated flow, it
%       can'be:
%         - 'd8' for the estimation of  the discrete flow directions given
%           by 8-connectivity, 
%         - 'gd8' for the estimation of the discrete global non-dispersive
%           flow directions proposed in [Paik08],
%         - 'ed82' and 'ed83' for the estimation of the discrete non-
%           dispersive flow directions of 2nd and 3rd orders also proposed 
%           in [Paik08];
%       default: flow='d8'.
%
% Property [propertyname  propertyvalues]
%   'dx', 'dy': optional horizontal and vertical pixel center spacing; if 
%       omitted, a value of 1.0 is assumed.
%   'mode' : see wflowacc.
%   'rflat' : see wflowacc.
%
% Outputs:
%   dir :  flow direction indices (sparse matrix).
%   slope : slope values (sparse matrix).
%   upslope : upslope area.
%
%
% Example from Paik:
%   dem = [ 93  99  96  95  94; ...
%           94  95  97  96  93 ; ...
%           91  98  100 97  95; ...
%           92  94  96  98  94; ...
%           89  91  90  92  93]
%
% References:
%   [Paik08]  K. Paik: "Global search algorithm for nondispersive flow path
%      extraction", Journal of Geophysical Research, 113:F04001, 2008.
%   [Tarb97]  D.G. Tarboton: "A new method for the determination of flow
%      directions and upslope areas in grid digital elevation models",
%      Water Resources Research, 33(2):309-319, 1997.
%   [WLD07]  J.P. Wilson, C.S. Lam and Y. Deng: "Comparison of the performance 
%      of flow-routing algorithms used in GIS-based hydrologic analysis",
%      Hydrological Processes, 21:1026-1044, 2007.

%% Parsing and checking parameters

error(nargchk(1, 14, nargin, 'struct'));
error(nargoutchk(1, 3, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(dem)
    error('grdsmooth:inputerror','a matrix is required in input'); 
end

p = createParser('SFLOW8');   % create an instance of the inputParser class.
% optional parameters
p.addOptional('flow', 'd8', @(x)ischar(x) && ...
    any(strcmpi(x,{'d8','gd8','ed82','ed83'})));
p.addParamValue('mode','default', @(x)ischar(x) && ...
    any(strcmpi(x,{'default','randomized','random'})));
p.addParamValue('rflat', false, @(x)islogical(x));
p.addParamValue('dx',1,@(x)isscalar(x) && round(x)==x && x>=1);
p.addParamValue('dy',1,@(x)isscalar(x) && round(x)==x && x>=1);
p.addParamValue('disp',false,@(x)islogical(x));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            


%% Checking variables and setting internal parameters

% general values
[NX,NY] = size(dem);
nrc = numel(dem);
% nans = isnan(dem);

[Y,X] = meshgrid(1:NY,1:NX);
% if any(size(X) ~= size(Y)) || any((size(X) ~= size(dem)));
%     error('X, Y and dem must have same size and same size as the dem')
% end


%% Computing the slopes

% find neighbours of cells: for each pixel in the dem, give the index list  
% of its 8 neighbours in the graph (unless the pixel is on any border of 
% the image)  
[ic1,icd1] = ixneighbours(dem);
% create the adjacency matrix 
AD = sparse(ic1,icd1,1,nrc,nrc);    

% calculate the local slope direction and maximum slope: for each pixel in the
% dem, give the slope between this pixel and its 8 neighbours (unless the pixel
% is on any border of the image)  
E = (dem(ic1)-dem(icd1)) ./ ...
    hypot((X(ic1)-X(icd1))/p.dx,(Y(ic1)-Y(icd1))/p.dy);
if nargout >= 2    % slope matrix
    varargout{1} = sparse(ic1,icd1,E,nrc,nrc);                     
end


%% Setting the (non unique) flow direction over the neighbours

% reset the null slopes
E(E<0) = 0;
% initialize the flow direction matrix with all slopes
M = sparse(ic1,icd1,E,nrc,nrc);

%% Routing through flats

% A flat or plateau exists when two or more adjacent cells have the same 
% elevation. Up to now flow direction indicates for these cells no exchange
% of water. The subsequent code first identifies flats and then iteratively
% removes them.

while p.rflat   
    % in a digital elevation model only those cells flats can be assigned as
    % flats that "do not give" and are -not- located on the edge of the dem.
    [r,ic,icd] = routeflat(dem,M,AD);
  
    if any(r)    
        % icd are the indices of the sills and ic are the indices of the
        % cells that contribute to the sills: water exchange from ic to icd
        ic  = ic(r);
        icd = icd(r);
    
        %  a new connects matrix is built for sills
        M = M+sparse(ic,icd,0.01,nrc,nrc);
        p.rflat = true;
    else
        p.rflat = false;
    end
end


%% Randomization of amount transferred to another cell

switch p.mode;
    case 'random'
        M = abs(sprandn(M));
    case 'randomized'
        % randomize coefficient. The higher, the more random
        rc = 0.01;
        M = M + (rc * abs(sprandn(M)));
    otherwise
end


%% Estimation of the flow direction given by direction of the maximal
% slope: LSD - local steepest slope

Ix = (1:nrc)';      
% retrieve the maxima values of the slopes and their position (index)
[maxS,Idx] = max(M,[],2); 
% get rid of the null values (if commented, set the variables i below when 
% defining Md and M2d)
ii = (maxS==0);
Ix(ii) = [];
Idx(ii) = [];
% M = sparse(Ix,Idx,1,nrc,nrc);

% dummy variables for checking that all pixels in the DEM have been visited
visited = zeros(size(Ix));

%% 'gd8','ed83' or 'ed82':  nondispersive single flow direction
if any(strcmp(p.flow,{'gd8','ed82','ed83'}))

    while ~all(visited(:))

        % local direction of the 1st LSD
        [Md,ix,iy,idx,idy] = localdownslope(Ix,Idx,NX,NY);
        % define the potential second direction pixels
        % possible choice for the 2nd direction: "the secondary direction is
        % defined as the steeper direction between clockwise and counterclockwise
        % directions adjacent to the LSD"
        A = localsecdir(ix,iy,idx,idy,NX,NY);
        [m, Idx2 ] = min( dem(A),[], 2 );        
        % final list of the 2nd LSD indices giving for each pixel (indexed
        % by Ix) the index of its secondary downslope direction
        Idx2 = diag(A(:,Idx2));
        % Idx2 : 2nd LSD indices giving for each pixel (indexed by Ix) the
        % index of its secondary downslope direction
        
        % extract all potential corridors (independently of the start cell):
        % "situation in which the flow passes through a
        % clearly defined corridor with steep lateral slopes: both cells in
        % clockwise and anticlockwise direciton adjacent to the LSD have
        % elevation higher than the LSD" (dem value on the 2nd downslope
        % direction > dem value on the LSD)"
        c = dem(Idx2) > dem(Ix);
        % get rid of those pixels as potential steepest slope and let the
        % flow pass in the corridor anyway (there wont be change for those
        % pixels)
        Idx2(c) = Idx(c);
        
       % select the start cell as the max of the dem
        [m,start] = max(dem(Ix).*(~visited));                        
        % and take the first maximal value in the dem we encounter
        start = start(1); % first found maximum
        % start : index position of the start cell in the vector Ix storing
        %   the position of the pixels with non null slope
        % define the seed pixel as the upstream cell 'directing' the flow
        upstream = start;
 
        % update the list of visited pixels
        visited(start) = 1;

        while ~isempty(start) && ~isempty(upstream)

            visited(upstream) = 1; 
             
            % look at the LSD: if a flat area is to be reached soon, then
            % stop propagating the flow for this start cell
            lsdupstream = find(Ix==Idx(upstream),1,'first');
            % lsdupstream : index position of the LSD of the upstream cell  
            % in the vector Ix 
            if isempty(lsdupstream) % we will reach a flat area
                upstream = [];                                         %#ok
                break; % we leave the LSD as it is
            end
            
            visited(lsdupstream) = 1;   
            
            if c(upstream) % then we are here in the presence of a corridor
                % these pixels will be used as new start cells for the
                % global search: "the LSD is chosen as the flow direction
                % and its downstream cell becomes a new start cell"
                start = lsdupstream;
                upstream = start;                                      %#ok
                % nothing else is changed, move downwards and go to the next
                % instance of the loop
                break;
           end % end of 'if c(upstream)...'
           % else proceed...
            
            % local direction of the 2nd LSD
            M2d = localdownslope(Ix,Idx2,NX,NY);
            
            % calculate the 'global' slope wrt the considered start cell
            % over the whole image
            % sstart = repmat(Ix(start),nrc,1);
            a = Ix(start);
            b = [Idx(lsdupstream) ;Idx2(lsdupstream)];
            SS = (dem(a)-dem(b)) ./ ...
                hypot((X(a)-X(b))/p.dx, (Y(a)-Y(b)/p.dy));
            
            % find where:
            %  - "the direction exists": Md(upstream)~=0 (always true),
            %  - "the flow is the same as the determined flow of the upstream
            %    cell": Md(lsdupstream)==Md(upstream),
            %  - "the secondary direction is the same as that of the upstream
            %    cell": M2d(lsdupstream)==M2d(upstream).
            % it gives the position of the set of pixels whose downslope 
            % direction could possibly be changed
            % refine by looking for the pixels where "the gradient between
            % the downstream cell of the secondary direction and the start
            % cell is steeper than the gradient between the downstream cell
            % of the primary direction and the start cell": SS(1)<SS(2)
            if Md(upstream) ~= 0 && ...               % existence
                    Md(lsdupstream) == Md(upstream) && ...  % unidirectionality
                    M2d(lsdupstream) == M2d(upstream) && ...% accrued deviation
                    SS(1) < SS(2) % gradient condition
               
                % update by replacing the final downstream direction by the
                % secondary direction
                Md(lsdupstream) = M2d(lsdupstream);
                Idx(lsdupstream) = Idx2(lsdupstream);
                % Idx2(ilsdupstream) = Idx2(Idx2(ilsdupstream));
                
            end % end of 'if Md(iupstream)...'
            
            % possibly update both the upstream and start cells by moving
            % them downwards, depending on the rule chosen for estimating
            % the flow
            switch p.flow
                case 'ed82'
                    % move the start cell downwards, 1 order apart from the
                    % cell of consideration
                    start = lsdupstream; % find(Ix==Idx(upstream),1,'first')
                    % consequently move the upstream cell
                    upstream = find(Ix==Idx(lsdupstream),1,'first'); % one
                    % position below because it has been already processed
                case  'ed83'
                    % move the start cell downwards, 2 orders apart from the
                    % cell of consideration
                    start = find(Ix==Idx(lsdupstream),1,'first');
                    upstream = start; % at that position, it has not been 
                    % processed
                otherwise % 'gd8'
                    % the start cell is left unchanged, the upstream cell
                    % is moved downwards
                    upstream = lsdupstream;
            end
            
        end % end of 'while ~isemtpy(iupstream)...'                      
    end % end of 'while sum(Ix)~=0...'

end

%% Sparse flow direction representation
% final sparse matrix giving, for each pixel, the index of the downslope 
% pixel in the matrix
M = sparse(Ix,Idx,1,nrc,nrc);

if p.disp
    displayflow(M,NX,NY);
end
    

%% Upslope area
if nargout ==3
    varargout{2} = upslopearea(dem,M);
end

end



%% ------------------------------------------------------------------------
function [r,ic,icd] = routeflat(dem,M,AD)

% check whether one of the surrounding cells is a giver. This is
% done by querying the neighbours of the flats.
ing = find(sum(M,2)==0);
ing(nans(ing)) = [];
a = full(sum(sparse(ing,ing,1,numel(dem),numel(dem))*AD,2)==8);
b = full(AD*a);

inb_flats = reshape(b<8 & a, size(dem));
IX_outb_flats = find(b & ~a);

% not too many of IX_outb_flats should be sills. To check whether a
% cell is a sill it
% 1. should be a giver and
% 2. have the same elevation as its neighbour in IX_inb_flats

[ic,icd] = ixneighbours(dem,inb_flats);
r   = (dem(ic)==dem(icd));
ic  = ic(r);
icd = icd(r);

r = ismembc(icd,IX_outb_flats);
end


%% ------------------------------------------------------------------------
function [Md,ix,iy,idx,idy] = localdownslope(Ix,Idx,NX,NY)

% column vector giving for each pixel the downslope direction in [1,9]
%  - for the pixel grid:
[ix,iy] = ind2sub([NX,NY],Ix);
%  - for the LSD directions:
[idx,idy] = ind2sub([NX,NY],Idx);

% find the local direction: the following mapping is used:
%          ---------------         -------------
%          | NW | N | NE |         | 1 | 4 | 7 |
%          ---------------         -------------
%          | W  |   | E  |    =>   | 2 |   | 8 |
%          ---------------         -------------
%          | SW | S | SE |         | 3 | 6 | 9 |
%          ---------------         -------------
% matrix giving for each pixel the downslope direction in [1,9]
Md = sub2ind([3,3],idx-ix+2,idy-iy+2);
% central pixel (5), no direction:   Md(Md==0)=5;
%  Md = zeros(nrc,1); % zeros(NX,NY);
%  %%in the case the variable ii above has not been set:
%  %%i = idx-ix+2>0 & idy-iy+2>0;
%  %%Md(Ix(i)) = sub2ind([3,3],idx(i)-ix(i)+2,idy(i)-iy(i)+2);
%  % Md = sub2ind([3,3],idx-ix+2,idy-iy+2);
%  Md(Ix) = sub2ind([3,3],idx-ix+2,idy-iy+2);
end

%% ------------------------------------------------------------------------
function neidir = localsecdir(ix,iy,idx,idy,NX,NY)
% coordinates of the clockwise neighbours of the LSD directions:
idx2 = idx + (idx<NX).*(ix>=idx).*(iy>idy) - ...
    (idx>1).*(ix<=idx).*(iy<idy);
idy2 = idy + (idy<NY).*(iy>=idy).*(ix<idx) - ...
    (idy>1).*(iy<=idy).*(ix>idx);
% column vector of corresponding indices:
neidir = sub2ind([NX,NY], idx2, idy2);
% note that Mclock is indexed by Ix

% coordinates of the counterclockwise neighbours of the LSD directions
idx2 = idx + (idx<NX).*(ix>=idx).*(iy<idy) - ...
    (idx>1).*(ix<=idx).*(iy>idy);
idy2 = idy + (idy<NY).*(iy>=idy).*(ix>idx) - ...
    (idy>1).*(iy<=idy).*(ix<idx);

% define the potential second direction pixels
neidir = cat(2, ...
    neidir, ...
    sub2ind([NX NY], idx2, idy2));

end


%% ------------------------------------------------------------------------
function Mflow = displayflow (M,x,y)
Mflow = zeros(x,y);
[uppix,downpix] = find(M);
[ix,iy] = ind2sub([x,y],uppix);
[idx,idy] = ind2sub([x,y],downpix);
Mflow(uppix) = sub2ind([3,3],idx-ix+2,idy-iy+2);

disp('flow coding:');
disp(reshape([1:4,0,6:9],3,3));

%disp('flow directions:');
%disp(Mflow);
figure, imagesc(Mflow), axis image, colormap jet

end


%% ------------------------------------------------------------------------
function mask = bordernans(dem)
% BORDERNANS - Find all connected NaNs connected to the DEM border.
%
%        mask = bordernans(dem);
%        
% Input: 
%    dem : input DEM image.
%
% Output:
%    mask : logical matrix the same size as dem; true values in mask 
%      correspond to border NaN locations. 
%
% Example:
%      E = magic(5);
%      E(2:5,1:2) = NaN;
%      E(3:4,4) = NaN
%      border_nans(E)

nan_mask = isnan(dem);
mask = nan_mask & ~imclearborder(nan_mask);
end


%% ------------------------------------------------------------------------
function A = upslopearea(dem,M)
% UPSLOPEAREA -  Compute the upslope area for each pixel of an input DEM
% given the flow matrix
%
%         A = upslopearea(dem,M);
%
% Inputs: 
%    dem : input DEM image.
%    M : (sparse) matrix representing the distribution of flow from pixel
%       to pixel.
%
% Output:
%    A : upslope area for each corresponding pixel of dem.
%
% Note: 
%  Connected groups of NaN pixels touching the border are treated as
%  having no contribution to flow.

% Right-side vector is normally all ones, reflecting an equal contribution
% to water flow originating in each pixel.
rhs = ones(numel(dem), 1);

% Connected groups of NaN pixels that touch the border do not contribute
% to water volume.
mask = bordernans(dem);
rhs(mask(:)) = 0;

A = M \ rhs;
A = reshape(A, size(dem));
end

