%% CONVOLUTION - Convolution based filtering.
%
%% Description
% Compute convolution with centered filter.
%
%% Syntax
%     y = CONVOLUTION(x, h);
%     y = CONVOLUTION(x, h, bound);
%
%% Inputs
% *|x|* : input signal; either 1D or 2D.
%
% *|h|* : centered filtered used to convolve |x|; it is centered at 0 for odd
%     length of the filter, and at 1/2 otherwise.
%
% *|bound|* : optional string parameter; it is either |'per'| for periodic 
%     extension of the signal or |'sym'| for symmetric extension.
%
%% Output
% *|y|* : convolved signal.
%
%% See also
% Related:
% <../../kernel/html/DIRGAUSSKERNEL.html |DIRGAUSSKERNEL|>,
% <../../kernel/html/GAUSSKERNEL.html |GAUSSKERNEL|>,
% <../../kernel/html/HOURGLASSKERNEL.html |HOURGLASSKERNEL|>,
% <matlab:webpub(whichpath('CONV')) |CONV|>,
% <matlab:webpub(whichpath('CONV2')) |CONV2|>,
% <matlab:webpub(whichpath('IFFT')) |IFFT|>,
% <matlab:webpub(whichpath('FFT')) |FFT|>,
% <matlab:webpub(whichpath('FFT2')) |FFT2|>.
% Called: 
% <CONVOLUTION_BASE.html |CONVOLUTION_BASE|>.

%% Function implementation
function y = convolution(x, h, varargin)

%%
% parsing parameters

error(nargchk(1, 11, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

if ~isnumeric(x)
    error('convolution:inputerror','numeric signal required in input');
end

p = createParser('CONVOLUTION');   % create an instance of the inputParser class.
p.addOptional('bound','sym', @(x)ischar(x) && any(strcmpi(x,{'sym','per'}))); 

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% main computation

[y,nd] = convolution_base(x, h, p.bound);

%%
% display

if p.disp
    figure,
    if nd==1,    plot(y,'r-');
    else
        imagesc(rescale(y,0,1)), axis image off;
        if size(y,3)==1,  colormap gray;  end
    end
    title('convolved signal');
end

end % end of convolution
