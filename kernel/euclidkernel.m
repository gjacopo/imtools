%% EUCLIDKERNEL - Euclidean distance weighted kernel.
%
%% Description
% Create a kernel spatially weighted with the (posibly inverse) Euclidean
% distance to a point (typically, belonging to the kernel) this is equivalent
% to a 'moving-window' style matrix with distance to centroid cell weighted.
% 
%% Syntax
%       W = EUCLIDKERNEL([wsx wsy]);
%       W = EUCLIDKERNEL([wsx wsy cind], [dy dx], norm, inverse );
%       W = EUCLIDKERNEL([wsx wsy cx cy], [dy dx], norm, inverse );
% 
%% Inputs
% *|wsx, wsy|* : vector of the (X,Y) dimension of the kernel window.
%
% *|cind|* or *|cx, cy|* : coordinates (passed as an integer >=1 matrix location
%     index or a couple of real >=0 (X,Y) position values) of the centroid 
%     of the kernel, ie. the point to which Euclidean distances are computed;
%     note that this point can also be set outside the kernel, however in the
%     case it is passed as an index, it is not garanteed that it will be
%     positioned as desired (no domain information is available); in the case
%     the coordinates of the centroid are not passed, it is set to a central
%     position of the kernel (with some arbitrary when the kernel dimensions 
%     are even).
%
% *|d|* : vector |[dy dx]| of the cellspacing in (X,Y) direction (typically,
%     |d=1|); default: |d=1|.
%
% *|norm|* : logical flag specifying if the weight matrix is to be normalized; 
%     default: |norm=false|.
%   
% *|inverse|* : logical flag setting if the inverse of the euclidean distance
%     is to be computed or not; default: |inverse=false|, ie. the weights 
%     are given by the Euclidean distance.
%
%% Output
% *|W|* : matrix with weights for every cell except center; weights are the 
%     (inverse of the) Euclidean distance to the center cell; center cell
%     weight is 0.
%
%% See also
% Related:
% <DIRGAUSSKERNEL.html |DIRGAUSSKERNEL|>,
% <GAUSSKERNEL.html |GAUSSKERNEL|>,
% <HOURGLASSKERNEL.html |HOURGLASSKERNEL|>.
% Called:
% <matlab:web(whichpath('MESHGRID')) |MESHGRID|>.

%% Function implementation
function W = euclidkernel(ws, d, n, inv)

%%
% parsing parameters

error(nargchk(1, 4, nargin, 'struct'));

% we allow for a variable number of entries
if nargin<4,  inv = false;
    if nargin<3,  n = false; 
        if nargin<2,  d = 1;  end
    end
end

%%
% setting parameters

if length(ws)<=2
    if length(ws)==1, ws = [ws ws];  end
    % if mod(ws(1),2)==0,  ws(1) = ws(1)+1;  end
    % if mod(ws(2),2)==0,  ws(2) = ws(2)+1;  end
    cx = ceil(ws(1)/2); cy = ceil(ws(2)/2);
    
elseif length(ws)==3
    [cx cy] = ind2sub([ws(1) ws(2)], ws(3));
    
elseif length(ws)==4
    cx = ws(3); cy = ws(4);
   
else 
    error('euclidkernel:inputerror',...
        'input dimension vector ws must be of dimension >=2 and <=4');
    
end

if cx>ws(1) || cy>ws(2)
    warning('euclidkernel:inputwarning',...
        'centroid outside kernel: its position may be different than what expected');
end

if length(d)==1
    d = [d d];   

elseif length(d)~=2
    error('euclidkernel:inputerror',...
        'input cellspacing vector must be of dimension <=2');
end

%%
% create a meshgrid centered around the passed centroid
[X,Y] = meshgrid(-cy+1:1:ws(2)-cy,-cx+1:1:ws(1)-cx);
X = X * d(2);
Y = Y * d(1);
W = sqrt(X.^2 + Y.^2);

%%
% inverse the weights
if inv
   W(cx,cy) = 1; %set center to 1 to avoid division by zero error
   W = 1./W; % get inverse distance
   W(cx,cy) = 0; %set center weight to zero
end

%% 
% normalize
if n
    s = sum(W(:));
    if s~=0,   W = W ./ s;  end
end

end % end of euclidkernel
