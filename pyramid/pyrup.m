%% PYRUP - Bottom-up hierarchical decomposition. 
% 
%% Description
% Perform a bottom-up hierarchical decomposition (resolution increases) of 
% an image into a pyramid using the function |EXPAND2D|; it also implements
% the classical bottom-up Laplacian Pyramid of [BA83].
% 
%% Syntax
%    Ip = PYRUP(I);
%    [Ip,Ep] = PYRUP(I, nlevels, stack, filter, order);
%
%% Inputs
% See function |PYRDOWN|.
%
% Note: default |nlevels=1|, ie. the pyramid has one level only.
% 
%% Outputs
% *|Ip|* : output pyramid; if |stack==true|, |Ip| is a cell of |nlevels| images, 
%     each being an expanded version of the previous level; if |stack=false|,
%     |Ip| is just the last level of the pyramid.
%
% *|Er|* : (optional) error of the last level of the pyramid, computed as
%     |error=signal­REDUCE2D(Ip)|. 
%
%% See also
% Related:
% <PYRDOWN.html |PYRDOWN|>, 
% <EXPAND2D.html |EXPAND2D|>, 
% <REDUCE2D.html |REDUCE2D|>.
% Called: 
% <PYRUP_BASE.html |PYRUP_BASE|>,
% <REDUCE2D_BASE.html |REDUCE2D_BASE|>.

%% Function implementation
function [Ip, varargout] = pyrup(I,varargin)

%%
% check if possible
error(nargchk(1, 13, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

if ~isnumeric(I)
    error('pyrup:inputparameter','a matrix is required in input'); 
end

%% 
% parsing parameters

p = createParser('PYRUP');   
p.addOptional('nlevels',1, @(x)isscalar(x) && x>=1);
p.addOptional('stack', true, @(x)islogical(x));
p.addOptional('filter','cspl', @(x)ischar(x) && ...
    any(strcmpi(x,{'lpl','spl','spl2','cspl','cspl2'}))); 
p.addOptional('order',3, @(x)isscalar(x) && x>0); 

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% checking parameters

if any(strcmpi(p.filter,{'cspl','cspl2''spl''spl2'}))
    if ~exist('expand2d_mex','file')
        error('pyrup:errorlibrary', ...
            'mex file expand2d_mex required for spline based pyramid');
    elseif isempty(p.order)
        p.order = 3;
    end
elseif ~isempty(p.order)
    warning('pyrup:warninginput', ...
        'input variable ''order'' ignored with Laplacian pyramid');
    p.order = [];
end

if (any(strcmpi(p.filter,{'cspl','cspl2'})) && ~ismember(p.order,[0,1,2,3,4])) || ...
        (strcmpi(p.filter,'spl') && ~ismember(p.order,[0,1,2,3])) || ...
        (strcmpi(p.filter,'spl2') && ~ismember(p.order,[0,1,3,5]))
    error('pyrup:inputerror', ...
        'check compatibility of filter and order variables');
end

%% 
% main computation

%%
% * compute the hierarchical representation
Ip = pyrup_base(I, p.nlevels, p.stack, p.filter, p.order);

%% 
% * compute the error
if nargout==2
    % take the last level
    if p.stack,   Io = Ip{p.nlevels}; 
    else          Io = Ip;       end;
    for i=1:p.nlevels
        [Io,Er] = reduce2d_base( Io, p.filter, p.order ); % expand
    end
    varargout{i} = Er;
end

end % end of pyrup