% -------------------------------------------------------------------------
function [ioF,varargout] = localorientfeature(F,varargin)
% INTERPORIENTFEATURE - 
%         ioF = localorientfeature(F);
%         [ioF, Zones, Omega] = localorientfeature(F, ...
%                  'Property', propertyvalue, ...);
%
% Inputs:
%   F : feature image with size [x, y, C] which values are locally
%      interpolated using the local (3x3) kernels defined by 'ker' or
%      'Kernel' and following the rules contained in 'Theta' or 'Omega' and
%      'Zones' (see below).
%
% Property [propertyname  propertyvalues]:
%   'ker' :
%   'Kernel' :
%   'Theta' :
%   'Omega' :
%   'Zones' :
%
% Outputs:
%   ioF : interpolated features.
%   Zones : optional image of zones computed when 'Theta' is passed.
%   Omega : optional interpolating factor computed when 'Theta' is passed.
%
% references:
%    [Leu00] J.G. Leu: "Edge sharpening through ramd width reduction",
%       Image and Vision Computing, 18: 501-514, 2000.

error(nargchk(1, 13, nargin, 'struct'));
error(nargoutchk(1, 3, nargout, 'struct'));

%% Prior checking and initialization
if ~isnumeric(F)
    error('localorientfeature:errorinput','a matrix is required in input');
end
[x, y, C] = size(F);
xy = x*y;

%% Parsing parameters
p = createParser('LOCALORIENTFEATURE');   
% only optional parameters
p.addParamValue('ker', [], @(x)ischar(x) && ...
    any(strcmpi(x,{'i0', 'ileu', 'i1', 'g0', 'gleu', 'g1'})));
p.addParamValue('Kernel', [], @(x)isnumeric(x) && size(x,1)==3 && size(x,2)==3);
p.addParamValue('Omega', [], @(x)isnumeric(x));
p.addParamValue('Zones', [], @(x)isnumeric(x));
p.addParamValue('Theta', [], @(x)isnumeric(x));
p.addParamValue('filt', 'mean', @(x)ischar(x) && ...
    any(strcmpi(x,{'mean','med'})));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            


%% Checking compatibility and setting default

% kernel and/or 'ker'
if isempty(p.Kernel) && isempty(p.ker)
    error('either ker or Kernel needs to be passe as an argument');
elseif ~isempty(p.ker) % ker should be a string
    if ~isempty(p.Kernel)
        disp('incompatible options Kernel and ker - Kernel ignored');
    end
    p.Kernel = local3x3kernel('ker',p.ker);
end
nl = size(p.Kernel,4);
nz = size(p.Kernel,3);

% orientation and zones
if ~isempty(p.Theta) && size(p.Theta,1)*size(p.Theta,2)~=xy
    error('input matrices Theta and F must have same (x-y) dimensions');  
elseif ~isempty(p.Omega) && size(p.Omega,1)*size(p.Omega,2)~=xy
    error('input matrices Omega and F must have same (x-y) dimensions');  
elseif ~isempty(p.Zones) && size(p.Zones,1)*size(p.Zones,2)~=xy
    error('input matrices Zones and F must have same (x-y) dimensions');
elseif (isempty(p.Omega) || isempty(p.Zones)) && isempty(p.Theta)
        error(['orientation Theta should be given when either Zones'...
        ' or Omega is not known']);  
end

%% Compute orientation and zones

if (isempty(p.Omega) || isempty(p.Zones))
    [p.Zones,p.Omega] = localorientzone(p.Theta,nz);
    if nargout>=2
        varargout{1} = p.Zones;
        if nargout==3, varargout{2} = p.Omega; end;
    end
end

%% Perform interpolation

% initialize the output
ioF = zeros(xy, nl, C);

% loop over C components
if strcmp(p.filt,'mean')
            p.filt=@imfilter;
elseif strcmp(p.filt,'med')
    p.filt = @wordfilt2;
end

for c=1:C
    ioF(:,:,c) = filtinterporient(p.filt, ...
         F(:,:,c), p.Zones, p.Omega, p.Kernel,nl,nz);
end

end
% end of localorientfeature


% -----------------------------------------------------------------------
function ioF = filtinterporient(filt,F,Zones,Omega,Kernel,nl,nz)
% initialize
ioF = zeros(length(Zones),nl);

% proceed for all zones and levels by filtering using the appropriate
% dependent (3 x 3) kernels stored in KernelPerZonePerLevel
for l = 1:nl
    % compute already the filtered version of the 1st zone, used as the
    % lower zone in the interpolation
     f = filt(F, Kernel(:,:,1,l));
    %f = imfilter(F, Kernel(:,:,1,l), 'replicate');
    for z = 1:nz
        % find the pixels whose orientation belongs to the current tested zone
        % (by keeping only the indexes of the pixels to the 'core matrix' (no
        % pixels on the border)
        ZZ = Zones(:) == z;
        % compute the filtered version for the next zone, corresponding to the
        % upper zone in the interpoation
         f1 = filt(F, Kernel(:,:,mod(z,nz)+1,l));
        %f1 = imfilter(F, Kernel(:,:,mod(z,nz)+1,l), 'replicate');
        % compute/update the (gradient or intensity) indices through
        % interpolation between adjacent orientations
        ioF(ZZ,l) = Omega(ZZ) .* (f(ZZ) - f1(ZZ)) + f1(ZZ);
        % note: there is a very low probability on natural images with
        % reasonable size that isempty(ZZ) is going to be true, therefore we
        % don't test it
        f = f1; % used in next loop
    end
end
end
% end of filtinterporient


% -----------------------------------------------------------------------
function ioF = spfiltinterporient(F,Zones,Omega,Kernel,nl,nz)

end
% end of spfiltinterporient
