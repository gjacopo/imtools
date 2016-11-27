%% UPSCALEXY - Image upsampling.
% 
%% Description
% Upsample monotically the input (possibly multispectral) image by either 
% duplicating, interpolating or inserting 0's between its input samples in 
% both X- and Y-directions.
% 
%% Syntax
%       upI = UPSCALEXY(I)
%       upI = UPSCALEXY(I, [sx sy], method)
%
%% Inputs
% *|I|* : 2D or 3D input matrix with size |(X,Y,Z)| to be upsampled in X-
%     (vertical) and Y- (horizontal) directions.
% 
% *|sx, sy|* : size (multiplicating factor) of the upsampling in X- and Y-
%     directions; default: |sx=sy=2|.
% 
% *|method|* : method used for upsampling, wich is either interpolation or
%     duplication; it can be either of:
% 
% * |'cubic'|, |'linear'|, |'nearest'|, |'splines'| or |'spline'| (call |INTERP2|)
%           if the matrix |I| is interpolated so that each block of the form
%           |[((i-1)*sx+1):(i*sx), ((j-1)*sy+1):(j*sy)]| of the output matrix 
%           is a quadratic interpolation between the values |I(i,j)|, |I(i+1,j)|,
%           |I(i,j+1)| and |I(i+1,j+1)| by the selected method.
% * |'replicate'| if the entries of the matrix are duplicated |sx| times
%           in X-direction and |sy| times in Y-direction so that each block
%           |[((i-1)*sx+1):(i*sx), ((j-1)*sy+1):(j*sy)]| of the output matrix
%           is filled with the identical value |I(i,j)|; creates a pixelized
%           image,
% * |'upsample'| to upsample the image |sx| times in X-direction and |sy| 
%            times in Y-direction by inserting 0 between the input samples,
%            similarly to the function |UPSAMPLE|;
% 
% default: |method='replicate'|.
%
%% Outputs
% *|upI|* : the upsampled matrix with size |(sx*X,sy*Y,Z)|.
%
%% See also 
% Related:
% <DOWNSCALEXY.html |DOWNSCALEXY|>, 
% <../../pyramid/html/EXPAND2D.html |EXPAND2D|>, 
% <matlab:webpub(whichpath('UPSAMPLE')) |UPSAMPLE|>,
% <matlab:webpub(whichpath('INTERP2')) |INTERP2|>.
% Called:
% <UPSCALEXY_BASE.html |UPSCALEXY_BASE|>, 

%% Function implementation
function upI = upscalexy(I,varargin)

%%
% parsing parameters

error(nargchk(1, 11, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

if ~isnumeric(I) && ~islogical(I)
    error('upscalexy:errorinput','matrix required in input');
end

% optional parameters
p = createParser('UPSCALEXY');   % create an instance of the inputParser class.
p.addOptional('S', [2 2 1], @(x)isnumeric(x) && all(x>=1) && length(x)<=3);
p.addOptional('method','replicate', @(x) ischar(x) && ... 
    any(strcmp(x,{'cubic','linear','spline','nearest','replicate','upsample'})));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            


%% 
% main computation: upscaling

upI = upscalexy_base(I, p.S, p.method);

%% 
% display

if p.disp
    figure, 
    [X,Y,C] = size(I);
    B = zeros(size(upI));
    for ic=1:C
        B(1:X,1:Y,ic) = I(:,:,ic);
    end
    subplot(1,2,1), imagesc(rescale(B,0,1)), axis image off, title('original');
    subplot(1,2,2), imagesc(rescale(upI,0,1)), axis image off, title('upscaled');
    if size(upI,3)==1,  colormap gray;  end
end

end % end of upscalexy


