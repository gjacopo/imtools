%% REDUCE2D - Reduce an image by a factor of 2 in each dimension.
% 
%% Description
% Wrapper for the |reduce2D| function of the [BIG] Spline Pyramids Software
% that reduces an image by a factor of 2 in each dimension using the
% algorithm described in [Unser99]; it also includes the code for reduction
% through the classical Laplacian Pyramid of [BA83].
%
%% Syntax
%    Ir = REDUCE2D(I);
%    [Ir, Er] = REDUCE2D(I, filter, order);
%
%% Inputs
% *|I|* : input image, possibly multispectral.
%
% *|filter|* : (optional) name of the filter; it is either:
%
% * |'lpl'| for Laplacian filter implemented in [BA83],
% * |'spl'| for the spline filter,
% * |'spl2'| for the spline $L^2$ filter,
% * |'cspl'| for the centered spline filter,
% * |'cspl2'| for the centered spline $L^2$ filter;
%
% for standard spline pyramids (|'spl'|, |'spl2'|), the coarser grid points
%       are  at the even integers; for centered pyramids (|'cspl'|, |'cspl2'|),
%       the coarser grid points are placed in-between their two closest finer
%       grid predecessors [BMIU99];  default: |filter='cspl'|.
%
% *|order|* : (optional) order for the filters based on splines; it is
%       either:
%
% * 0, 1, 2 or 3 for the spline filter,
% * 0, 1, 3 or 5 for the spline L2 filter,
% * 0, 1, 2, 3 or 4 for the both centered spline filters,
%
% in the case |filter='lpl'|, this parameter is ignored; default: |order=3|,
% with |filter='cspl'|.
%
%% Outputs
% *|Ir|* : image reduced by a factor 2.
%
% *|Er|* : (optional) error image, computed as |error=signal­EXPAND2D(Ir);|
%       it gives the loss of information due to image reduction. 
%
%% Acknowledgment
% This software package implements the basic |REDUCE| and |EXPAND| operators 
% for the reduction and enlargement of signals and images by factors of two.
% A signal is represented by a polynomial spline which is a continuously-
% defined function [Unser99]; the model is interpolating and is unambiguously 
% defined by the sample values. The spline model specifies the enlargement
% mechanism (spline interpolation), as well as the reduction algorithm which 
% is optimal in the least squares sense. 
%
%% References
% [BA83]  P.J. Burt and E.H. Adelson: "The Laplacian Pyramid as a Compact
%      Image Code", IEEE Transactions on Commununication, 31(4):337-345, 1983.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1095851&tag=1>
%
% [Unser99]  M. Unser: "Splines: a perfect fit for signal and image 
%      processing", IEEE Signal Processing Magazine, 16(6):22-38, 1999.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=799930>
%
% [BMIU99]  P. Brigger, F. Muller, K. Illgner and M. Unser: "Centered
%      pyramids," IEEE Transactions on Image Processing, 8(9):1254-1264, 
%      1999.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=784437>
%
% [BIG]  Software available at:      http://bigwww.epfl.ch/sage/pyramids/.
%      See README at:        http://bigwww.epfl.ch/sage/pyramids/README.TXT
%
%% See also
% Related:
% <PYRDOWN.html |PYRDOWN|>, 
% <EXPAND2D.html |EXPAND2D|>, 
% <../../statistics/html/DOWNSCALEXY.html |DOWNSCALEXY|>.
% Called: 
% <REDUCE2D_BASE.html |REDUCE2D_BASE|>, 
% <EXPAND2D_BASE.html |EXPAND2D_BASE|>.

%% Function implementation
function [Ir,varargout] = reduce2d(I, varargin)

%%
% check if possible
error(nargchk(1, 11, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

if ~isnumeric(I)
    error('reduce2d:inputparameter','a matrix is required in input'); 
end

%% 
% parsing parameters

p = createParser('REDUCE2D');   
p.addOptional('filter','cspl', @(x)ischar(x) && ...
    any(strcmpi(x,{'spl','spl2','cspl','cspl2'}))); 
p.addOptional('order',[], @(x)isempty(x) || (isscalar(x) && x>0 && x<=5)); 

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% checking parameters

if any(strcmpi(p.filter,{'cspl','cspl2''spl''spl2'}))
    if ~exist('reduce2d_mex','file')
        error('reduce2d:errorlibrary', ...
            'mex file reduce2d_mex required for spline based pyramid');
    elseif isempty(p.order)
        p.order = 3;
    end
elseif ~isempty(p.order)
    warning('reduce2d:warninginput', ...
        'input variable ''order'' ignored with Laplacian pyramid');
    p.order = [];
end

if (any(strcmpi(p.filter,{'cspl','cspl2'})) && ~ismember(p.order,[0,1,2,3,4])) || ...
        (strcmpi(p.filter,'spl') && ~ismember(p.order,[0,1,2,3])) || ...
        (strcmpi(p.filter,'spl2') && ~ismember(p.order,[0,1,3,5]))
    error('reduce2d:inputerror', ...
        'check compatibility of filter and order variables');
end

%%
% main calculation

Ir = reduce2d_base(I, p.filter, p.order); 

if nargout==2
    J = expand2d_base(Ir, p.filter, p.order);
    varargout{1} = abs(I - J);
end

%%
% display

if p.disp
    figure, subplot(1,nargout,1)
    imagesc(rescale(Ir,0,1)), axis image off, title('reduced image');
    if size(I,3)==1, colormap gray;  end;
    if nargout==2
        subplot(1,2,2)
        imagesc(rescale(varargout{1},0,1)), axis image off, title('error image');
        if size(I,3)==1, colormap gray;  end;
    end
end

end % end of reduced2d