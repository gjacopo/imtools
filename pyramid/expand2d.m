%% EXPAND2D - Upsample a signal by a factor of 2.
% 
%% Description
% Wrapper for the |expand2d| function of the [BIG] Spline Pyramids Software
% that upsamples a signal by a factor of 2 and then fills-in the missing
% samples by application of a suitable interpolation filter. 
% The interpolation filter is specified by the underlying spline model.
% It also implements the code for expansion through the classical Laplacian
% Pyramid of [BA83]
%
%% Syntax
%    Ie = EXPAND2D(I);
%    [Ie, Ee] = EXPAND2D(I, filter, order);
%
%% Inputs
% See function |REDUCE2D|.
%
%% Outputs
% *|Ir|* : image expanded by a factor 2.
%
% *|Er|* : (optional) error image, computed as |error=signal­REDUCE2D(Ir)|. 
% 
%% See also
% Related:
% <PYRUP.html |PYRUP|>, 
% <REDUCE2D.html |REDUCE2D|>, 
% <../../geometry/html/UPSCALEXY.html |UPSCALEXY|>.
% Called: 
% <REDUCE2D_BASE.html |REDUCE2D_BASE|>, 
% <EXPAND2D_BASE.html |EXPAND2D_BASE|>.

%% Function implementation
function [Ie,varargout] = expand2d(I, varargin)

%%
% check if possible
error(nargchk(1, 3, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

if ~isnumeric(I)
    error('expand2d:inputparameter','matrix required in input'); 
end


%% 
% parsing parameters

p = createParser('EXPAND2D');   
p.addOptional('filter','cspl', @(x)ischar(x) && ...
    any(strcmpi(x,{'spl','spl2','cspl','cspl2'}))); 
p.addOptional('order',3, @(x)isempty(x) || (isscalar(x) && x>0 && x<=5)); 

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            


%% 
% checking parameters
if any(strcmpi(p.filter,{'cspl','cspl2''spl''spl2'}))
    if ~exist('expand2d_mex','file')
        error('expand2d:errorlibrary', ...
            'mex file expand2d_mex required for spline based pyramid');
    elseif isempty(p.order)
        p.order = 3;
    end
elseif ~isempty(p.order)
    warning('expand2d:warninginput', ...
        'input variable ''order'' ignored with Laplacian pyramid');
    p.order = [];
end

if (any(strcmpi(p.filter,{'cspl','cspl2'})) && ~ismember(p.order,[0,1,2,3,4])) || ...
        (strcmpi(p.filter,'spl') && ~ismember(p.order,[0,1,2,3])) || ...
        (strcmpi(p.filter,'spl2') && ~ismember(p.order,[0,1,3,5]))
    error('expand2d:inputerror', ...
        'check compatibility of filter and order variables');
end


%% 
% main calculation: dealing with multispectral images

Ie = expand2d_base(I, p.filter, p.order); 

if nargout==2
    J = reduce2d_base(Ie, p.filter, p.order);
    varargout{1} = abs(I - J);
end

%%
% displaying 

if p.disp
    figure, subplot(1,nargout,1)
    imagesc(rescale(Ie,0,1)), axis image off, title('expanded image');
    if size(I,3)==1, colormap gray;  end;
    if nargout==2
        subplot(1,2,2)
        imagesc(rescale(varargout{1},0,1)), axis image off, title('error image');
        if size(I,3)==1, colormap gray;  end;
    end
end

end % end of expand2d