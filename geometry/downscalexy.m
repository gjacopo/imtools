%% DOWNSCALEXY - Image downsampling. NOT IMPLEMENTED YET.
% 
%% Description
% Downsample monotically the input matrix by keeping every n-th sample in
% both directions starting with the first. 

%% Syntax
%     downI = DOWNSCALEXY(I);
%     downI = DOWNSCALEXY(I, [sx sy]);
%
%% Inputs
% *|I|* : 2D or 3D input matrix with size |(X,Y,Z)| to be downsampled in X-
%     (vertical) and Y- (horizontal) directions.
% 
% *|sx, sy|* : size (dividing factor) of the downsampling in X- and Y-
%     directions; default: |sx=sy=2|.
% 
% *|sy|* : ibid in ; default: |sy=sx|.
%
%% Outputs
% *|downI|* : the downsampled matrix with size |(X/sx,Y/sy,Z)|.
%
%% Credit
% <mailto:grazzja@lanl.gov J.Grazzini> (ISR-2/LANL)
%
%% See also 
% Related:
% <UPSCALEXY.html |UPSCALEXY|>, 
% <matlab:webpub(whichpath('DOWNSAMPLE')) |DOWNSAMPLE|>.

%% Function implementation
function downI = downscalexy(I,varargin)                               %#ok                        

error('downscalexy:methoderror', 'method not yet implemented');

%%
% parsing parameters
error(nargchk(1, 2, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

if ~isnumeric(I)
    error('downscalexy:errorinput','a matrix is required in input');
end

% optional parameters
p = createParser('DOWNSCALEXY');   % create an instance of the inputParser class.
p.addOptional('S', [2 2 1], @(x)isnumeric(x) && all(x>=1) && length(x)<=2);
% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p); 

%% 
% internal variables

%% 
% downscaling

end % end of downscalexy


