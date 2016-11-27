%% EXTREMA3X3 - Basic 2D statistical operations in local 3x3 neighbourhoods.
%
%% Description
% Perform 2D erosion/dilation operations of an image in 3x3 neighbourhoods
% (with connectivity 4 or 8) using standard 1D extrema (min/max) operations.
% Useful when the Image Processing Toolbox is not available.
%
%% Syntax
%        M = extrema3x3(I, fext, nconn);
%
%% Inputs
% *|I|* : input image, possibly multispectral.
%
% *|fext|* : function handle; it is either @min (for erosion) or @max (for
%     dilation); default: |fext=@min|, hence the image is eroded.
%
% *|nconn|* : connectivity index; it is either 4 or 8; default: |conn=8|.
%
%% Outputs
% *|M|* : eroded or dilated image with same dimension as input image |I|.
%
%% See also
% Related:
% <matlab:webpub(whichpath('IMERODE')) |IMERODE|>,
% <matlab:webpub(whichpath('IMDILATE')) |IMDILATE|>.

%% Function implementation
function I = extrema3x3(I, varargin)

%% 
% parsing parameters
narginchk(1, 3);
nargoutchk(1, 1);
if ~isnumeric(I)
    error('extrema3x3:inputerror','a matrix is required in input');
end

p = createParser('EXTREMA3X3');   % create an instance of the inputParser class.
p.addOptional('fext', @min, @(x) isa(x,'function_handle') && ...
    (isequal(x,@min) || isequal(x,@max)));
p.addOptional('nconn', 8, @(x) (x==4 || x==8));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%if isequal(p.fext,@min) || isequal(p.fext,@max)
p.fext = @(x,y) p.fext(x,[],y);
%else
%p.fext = @(x,y) p.fext(x,y);
%end

%%
% distinguish special multispectral case
C = size(I,3);

if C>1
    ext = zeros(size(I)); 
    for c=1:C
        ext(:,:,c) = extrema3x3(I(:,:,c), p.fext, p.nconn);
    end
    return;
end

%%
% define 1D (columnwise) basic separable operation

%   %----------------------------------------------------------------------
    function [I,m,n] = colextrema3x3(I,fext)
        [m,n] = size(I);
        I = padarray(I,[0,1],'symmetric','both');
        I = reshape(fext(...
            cat(2, reshape(I(:,2:end-1),m*n,1), ...
                reshape(I(:,1:end-2),m*n,1), ...
                reshape(I(:,3:end),m*n,1)), ...
                2),m,n);
    end
%   %----------------------------------------------------------------------

%% 
% perform local neighbourhood operation
[A, m, n] = colextrema3x3(I,p.fext);
if p.nconn==8
    I = colextrema3x3(A',p.fext)';
elseif p.nconn==4
    I = colextrema3x3(I',p.fext)';
    I = reshape(p.fext(cat(2,A(:),I(:)),2),m,n);
end
    
end
