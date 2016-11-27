%% BLURMAP - Adaptive local isotropic Gaussian bluring.
%
%% Description
% Perform an adaptive local isotropic Gaussian bluring of an image where the
% scale of the Gaussian (standard deviation) is given by a blur map like in
% [EZ98].
%
%% Syntax
%     Iblur = BLURMAP(I);
%     [Iblur, qmap] = BLURMAP(I, blur, qstep);
%
%% Inputs
% *|I|* : an input image with size |(X,Y,C)|, where |C>1| when |I| is multichannel
%     (in that latter case, the same bluring is applied over the different 
%     channels separately).
%
% *|blur|* : optional image setting locally (on each pixel) the bluring 
%     parameter to be used in the isotropic Gaussian smoothing; it must be 
%     either: 
%
% * a scalar, then uniform isotropic bluring is performed all over the 
%         input image, 
% * a row (resp. colum) vector of size |X| (resp. size |Y|), then a row
%         (resp. column) dependent bluring with sigma parameters constant
%         over the columns (resp. rows) is performed, similarly to the test
%         image used in [EZ98]), 
% * a matrix of same size |(X,Y)| as the input matrix, where each entry
%         indicates the standard deviation of the Gaussian filter to be 
%         applied at each pixel;
%
% default: the blur values are normally distributed with mean 2 and standard 
%     deviation 0.5; note moreover that the map is quantized to avoid a too
%     large number of possible sigma values.
%
% *|qstep|* : optional quantizing parameter.
% 
%% Outputs
% *|Iblur|* : adaptively locally blurred image.
%
% *|qmap|* : optional output image storing the quantized local sigma values
%     actually used in the bluring.
%
%% Example
%   a=zeros(256,256);
%   a(:,128:256)=255;
%   f=linspace(0.25,10,size(a,1))';
%   s=f(2)-f(1);
%   b=blurmap(a,f,s); % bluring similarly to [EZ98]
%   figure, imagesc(b), colormap gray
% 
%% Reference
% [EZ98] J.H. Elder and S.W.Zucker: "Local scale control for edge
%      detection and blur estimation", IEEE Trans. on Pattern Analysis and
%      Machine Intelligence, 20:699-716, 1998.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=689301&tag=1>
%
%% See also
% Related:
% <CONVOLUTION.html |CONVOLUTION|>,
% <ADAPTIVEFILT.html |ADAPTIVEFILT|>.
% Called: 
% <BLURMAP_BASE.html |BLURMAP_BASE|>.

%% Function implementation
function [Iblur,varargout] = blurmap(I, varargin)

%%
% parsing parameters

% mandatory parameter
if ~isnumeric(I)
    error('blurmap:inputerror','a matrix is required in input');
end

% check number of parameters
error(nargoutchk(1, 2, nargout, 'struct'))
error(nargchk(1, 3, nargin, 'struct'))

% optional parameters
p = createParser('BLURMAP');   
% principal optional parameters
p.addOptional('blur', [], @(x)isnumeric(x) && ~any(x(:))<0.05);
p.addOptional('qstep', .3, @(x)isscalar(x) && isfloat(x) && x>0);

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%% 
% checking and setting parameters

[X,Y] = size(I(:,:,1));

% compute the bluring map (note that by doing it now, in the case of
% multispectral image, the same bluring map will be used for buring all the 
% channels.
if isempty(p.blur)
    p.blur = 2 + randn(X,Y)/2;
else
    [x,y] = size(p.blur); % other dimensions are ignored
    s = min(x,y); S = max(x,y);
    if s == 1
        if (s==x && S~=Y) ||(s==y && S~=X)
            error('blurmap:inputerror',...
                'input vector must be consistent with input image');
        else
            p.blur = repmat(p.blur,(s==x)*X+(s~=x),(s==y)*Y+(s~=y));
        end
    elseif x ~= X || y ~= Y
        error('blurmap:inputerror','input matrices must have same size');
    end
end

% dealing with multispectral images
% if C>1
%     Iblur = zeros(X, Y, C);
%     if nargout==2,  varargout{1} = zeros(X, Y, C);    end
%     for c=1:C
%         [Iblur(:,:,c) tmp] = blurmap(I(:,:,c), 'parse', p);
%     end
%     if nargout==2,  varargout{1}(:,:,c) = tmp;    end     
%     return;
% end

%% 
% main computation 

[Iblur,quantblur] = blurmap_base(I, p.blur, p.qstep);

if nargout==2
    varargout{1} = quantblur;
end

%%
% display

if p.disp
    figure, imagesc(rescale(Iblur,0,1)), axis image off, title('blurred image');
end

end % end of blurmap
