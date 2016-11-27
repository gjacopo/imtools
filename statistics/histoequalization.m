%% HISTOEQUALIZATION - Histogram equalization of an image.
% 
%% Description
% Perform histogram equalization, ie. change the values of the input image
% so that its reordering match the ordered values of another vector/matrix 
% or a given density.
%
%% Syntax
%   X = HISTOEQUALIZATION(x, y);
%   X = HISTOEQUALIZATION(x, y, 'Property', propertyvalue, ...);
%   X = HISTOEQUALIZATION(x, 'gauss'/'lin', 'Property', propertyvalue, ...);
%
%% Inputs
% *|x|* : vector/matrix to equalize.
%
% *|y|* : reference vector/matrix to operate the equalization of x, or string
%       set to |'lin'| or |'gauss'| for linear or gaussian equalization resp.
%
%% Property [propertyname  propertyvalues]
% *|'dir'|* : string set to operate the equalization along a preferential 
%       direction when the input is a multidimensional matrix; it can be
%       either |'col'|, |'row'| or |'vec'| to consider the distribution along
%       the columns, the rows or the 3rd dimension resp.; default: |dir=[]|,
%       ie. no preferential direction for equalizing the values, the
%       equalization is performed wrt the distribution of all the values in
%       the matrix |x|.
%
% *|'abs'|* : logical flag set to |true| to operate only on absolute values;
%       default: |abs=false|.
%
% *|'lam'|* : scalar in |[0,1]| to interpolate between the two histograms;
%       default: |lam=1|.
%
%% Outputs
% *|X|* : equalized version of |x|.
%
%% See also
% Related:
% <matlab:webpub(whichpath('HIST')) |HIST|>.
% Called:
% <HISTOEQUALIZATION_BASE.html |HISTOEQUALIZATION_BASE|>.

%% Function implementation
function X = histoequalization(x, varargin)

%%
% parsing parameters

narginchk(1, 16);
nargoutchk(1, 1);

if ~isnumeric(x)
    error('histoequalization:inputparameter','vector or matrix required in input'); 
end

p = createParser('HISTOEQUALIZATION');   
p.addRequired('y', @(x)isnumeric(x) || ...
    (ischar(x) && any(strcmpi(x,{'lin','gauss'}))));  % has to be entered!!!
p.addParamValue('dir',[], @(x)isempty(x) || ...
    (ischar(x) && any(strcmpi(x,{'row','col','vec'})))); 
p.addParamValue('lam', 1, @(x)isscalar(x) && x<=1 && x>=0);
p.addParamValue('abs', false, @(x)islogical(x));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% checking variables

if ~isempty(p.dir) && ischar(p.dir) && isnumeric(p.y)
    if strcmpi(p.dir,'col') && ~size(x,2)>1 || ~(size(x,2)==size(p.y,2))
        error('histoequalization:inputerror', ...
            'x and y must have same number of columns');
        
        
    elseif strcmpi(p.dir,'row') && ~size(x,1)>1 || ~(size(x,1)==size(p.y,1))
        error('histoequalization:inputerror', ...
            'x and y must have same number of rows');
        
    elseif strcmpi(p.dir,'vec') && ~size(x,3)>1 || ~(size(x,3)==size(p.y,3))
        error('histoequalization:inputerror', ...
            'x and y must have same number of dim 3');
    end
end

%%
% main computation

X = histoequalization_base(x, p.y, p.dir, p.abs, p.lam);

end % end of histoequalization
