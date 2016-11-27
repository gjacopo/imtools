%% FINDLOCALEXTREMA - Local neighbourhood extrema of a map.
%
%% Description
% Extract the extrema (min or max) of any multispectral image in local square
% neigbourhoods and outputs a logical map.
%
%% Syntax
%     M = FINDLOCALEXTREMA(I);
%     M = FINDLOCALEXTREMA(I, ord, win, 'Property', propertyvalue, ... );
%
%% Inputs
% *|I|* : input image of size [X,Y,C], possibly multichannel when C>1.
%   rank : string variable defining the extrema to be computed locally; it 
%     is either 'min' or 'max' to compute local minima or maxima resp.
%   
% *|win|* : width of the window of analysis.
%
%% Property [propertyname  propertyvalues]
% *|'graph'|* : parameter defining the implicit graph connectivity assumed when
%     computing the local rank filter; it is either 4 or 8; default: |graph=8|.
%
% *|'corr'|* : string defining the way the channels of the input image are
%     combined in the output extrema map when |C>1|; it can be:
%
% * |'forall'| if the pixels to be output are the extrema considered
%          separately in the different channels,
% * |'inall'| if the pixels are the extrema simultaneously found in all
%          channels
% * |'inone'| if the pixels are extrema in one channel at least;
%     
% |corr| is naturally ignored with grayscale images; default: |corr='forall'|.
%
%% Outputs
% *|M|* : logical output map where extrema are assigned true; depending on the
%     input variable |corr| (see above), |M| can be of size (|X,Y|) (cases
%     |'inall'| and |'inone'|) or of size |(X,Y,C)| (case |'forall'|).
% 
%% See also
% Related:
% <FINDLOCALMAX.html |FINDLOCALMAX|>,
% <../../sharpen/html/MAPTRANSITION.html |MAPTRANSITION|>,
% <matlab:webpub(whichpath('ORDFILT2')) |ORDFILT2|>.
% Called: 
% <FINDLOCALEXTREMA_BASE.html |FINDLOCALEXTREMA_BASE|>.

%% Function implementation
function M = findlocalextrema(I,varargin)

%%
% parsing parameters
error(nargchk(1, 15, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(I)
    error('findlocalextrema:inputerror','a matrix is required in input'); 
end

% optional parameters
p = createParser('FINDLOCALEXTREMA');   % create an instance of the inputParser class.
p.addOptional('ord', 'max', @(x) ischar(x) && any(strcmpi(x,{'min','max'})));
p.addOptional('win',3, @(x)isscalar(x) && isfloat(x) && x>=3);
p.addParamValue('graph',8, @(x)isscalar(x) && (x==4 || x==8));
p.addParamValue('corr',[], @(x)ischar(x)&& ...
    any(strcmpi(x,{'forall','inall','inone'})));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%% 
% setting internal variables

C = size(I,3);

if isempty(p.corr),   p.corr = 'forall';
elseif C==1 
    warning('findlocalextrema:inputwarning', 'variable corr ignored with scalar image');
end

%% 
% main computation

M = findlocalextrema_base(I, p.ord, p.win, p.graph, p.corr);

%%
% display

if p.disp
    figure, imagesc(rescale(M,0,1)), axis image off, title('local extrema')
    if ~strcmp(p.corr,'forall'), colormap gray; end;
end

end % end of findlocalextrema
