function flowacc = sfdacc(dem,M)

% SFDACC - Upslope area algorithm for Digital Elevation Models using
% single flow direction
%
%      flowacc = sflowacc(dem,M)
%      flowacc = sflowacc(dem,M,'propertyname',propertyvalue,...)
% 
% Single flow accumulation algorithm that routes through flat terrain (not 
% sinks). Remove sinks with the function imfill that requires the image 
% processing toolbox.
%
% Input: 
%   dem       digital elevation model same size as X and Y
%   M         flow direction computed using any single flow direction
%             technique
%
% Properties:
% propertyname               propertyvalues
%   'W0'              see wflowacc
%   'edges'           see wflowacc
%
% Output:
% flowacc   flowaccumulation (upslope area) grid


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

% check input using PARSEARGS
params.W0          = ones(siz);
params.edges       = {'closed','open'};
params = parseargs(params,varargin{:});

% ********************************************************************

% find neighbours of cells: for each pixel in the dem, give the index list
% of its 8 neighbours in the graph (unless the pixel is on any border of 
% the image)  
[ic1,icd1] = ixneighbours(dem);
% create the adjacency matrix 
AD = sparse(ic1,icd1,1,nrc,nrc);    

% edge correction when open
switch params.edges
    case 'open';
        edgecorr = full(sum(AD,2)/8);
end

% ******************************************************************
% defining the accumulation 

switch params.edges
    case 'open';
        flowacc = (speye(nrc,nrc)-spdiags(edgecorr,0,nrc,nrc)*M')\params.W0(:);
    otherwise
        flowacc = (speye(nrc,nrc)-M')\params.W0(:);        
end

flowacc = reshape(flowacc,siz);

% and this is it...
return;



