%% ASF - Alternate Sequential Filtering.
%
%% Description
% Perform the Alternate Sequential Filtering (ASF) and the ASF by reconstruction
% of an image [Serra88,Soille03].
%
%% Syntax
%     F = ASF(I);
%     F = ASF(I, op, n);
%     F = ASF(I, op, n, shape);
%     F = ASF(I, op, n, shape, s1);
%     F = ASF(I, op, n, shape, s1, s2);
%
%% Inputs
% *|I|* : input image, possibly multispectral.
%
% *|op|* : string defining if the filter is by reconstruction (it should start
%     with a |'r'| then), and specifying the order in which the opening (|'o'|)
%     and closing (|'c'|) operators are applied; it can be either:
%    
% * |'oc'| or |'co'| for standard ASF applying resp. first the opening, then
%         the closing, or first the opening, then the closing,
% * |'oco'| or |'coc'| for standard ASF applying resp., in this order, 
%         opening, closing and opening or closing, opening and closing,
% * |'rco'| or |'roc'| for ASF by reconstruction applying resp. first the
%         opening, then the closing, or first the closing, then the opening, 
% * |'roco'| or |'rcoc'| for ASF by reconstruction applying resp., in this
%         order, opening, closing and opening or closing, opening and 
%         closing, 
%    
% default: |op='oc'|.
%
% *|n|* : number of iterations of the sequential filter; default: |n=1|.
%
% *|shape|* : any string among |'disk'|, |'rectangle'|, |'square'|, |'diamond'|, |'line'|,
%     |'periodicline'| and |'octagon'| (in particular, |'ball'|, 'arbitrary' and 
%     |'pair'| are not supported) or any scalar value. In the former case, a 
%     default SE of type |shape| and size 3 will be output, in the latter 
%     case, a SE of type |'disk'| and size shape will be output; default: 
%     |shape='disk'|.
%
% *|s1|* : optional shape (string, see above) or size of the SE when |nargin=2|;
%     when |nargin=3| (ie. |s2| is also defined), |s1| is necessarly a scalar
%     value.
%
% *|s2|* : optional scalar value for defining the SE (see |STREL| possibly defined
%     with 2 parameters: |'periodicline'|, |'disk'|, |'line'|).
%
%% Outputs
% *|F|* : output of the ASF or the ASF by reconstruction applied on the input
%     image. 
%
%% References
% [Serra88]  J. Serra: "Image Analysis and Mathematical Morphology", vol. 2, 
%      Academic Press, New York (USA), 1988. 
%
% [Soille03]  P. Soille: "Morphological Image Analysis?Principles and 
%      Applications", 2nd ed. Berlin (Germany), Springer-Verlag, 2003. 
%
%% See also
% Related:
% <matlab:web(whichpath('IMOPEN')) |IMOPEN|>,
% <matlab:web(whichpath('IMCLOSE')) |IMCLOSE|>,
% <matlab:web(whichpath('IMRECONSTRUCT.html')) |IMRECONSTRUCT|>,
% <IMROPEN.html |IMROPEN|>,
% <IMRCLOSE.html |IMRCLOSE|>,
% <FLATSTREL.html |FLATSTREL|>.
% Called: 
% <ASF_BASE.html |ASF_BASE|>.

%% Function implementation
function F = asf(I, varargin)

if isempty(ver('images'))
    error('asf:errortoolbox', 'Image Processing toolbox required');
end

error(nargchk(1, 17, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

%%
% parsing parameters

if ~isnumeric(I)
    error('asf:inputparameter','a matrix is required in input'); 
end

% optional parameters
p = createParser('ASF');   

p.addOptional('op', 'oc', @(x)ischar(x) && ...
    any(strcmpi(x,{'co','oc','oco','coc','rco','roc','rcoc','roco'})));
p.addOptional('n', 1, @(x)isnumeric(x) && x>=1);
p.addOptional('shape', 'disk', @(x) ischar(x) && ...
    any(strcmpi(x,{'disk','rectangle','square','diamond', ...
                'line','periodicline','octagon'})));
p.addParamValue('s1', 3, @(x)isnumeric(x));
p.addParamValue('s2', [], @(x)isnumeric(x));

% parse and validate all input arguments
p.parse(varargin{:});
p = getvarParser(p);

%% 
% main processing

F = asf_base(I, p.n, p.shape, p.op, p.s1, p.s2);

%% 
% display

if p.disp
    thetitle = 'ASF';
    if p.op(1)=='r', thetitle = [thetitle ' by reconstruction'];  end
    thetitle = [thetitle ' - order: ' p.op(end-1:end)];
    figure, imagesc(rescale(F)), axis image off, title(thetitle);
    if size(I,3)==1,  colormap gray;  end
end

end % end of asf
