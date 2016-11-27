%% GRANULOMETRY - Morphological granulometry representation of an image.
%
%% Description
% Compute a granulometry through a series of morphological operations (filters)
% following the various approaches described in [Mathe75,Serra82,CD94,DC99].
%
%% Syntax
%     G = GRANULOMETRY(I);
%     G = GRANULOMETRY(I, op);
%     G = GRANULOMETRY(I, op, 'Property', propertyvalue, ... );
%
%% Inputs
% *|I|* : input image of size |(X,Y,C)|, with |C>1| when |I| is multispectral.
%
% *|op|* : (optional) string setting the elementary operations performed at
%     each step for computing the final granulometry; it can any of those
%     strings: |'o'| (or |'open'|), |'c'| (or |'close'|), |'oc'|, |'co'|,
%     |'oco'|, |'coc'|, |'ro'|, |'rc'|, |'roc'|, |'rco'|, |'rcoc'|, |'roco'|,
%     |'e'| (or |'erode'|) or |'d'| (or |'dilate'|), where the |'o'| represents
%     the opening operation (see function |IMOPEN|), the 'c' represents the
%     closing operation (see function |IMCLOSE|), the prior 'r' means that
%     the reconstruction version of these operators are considered (see
%     functions |IMROPEN| and |IMRCLOSE|), and where the order of the letters
%     represents the order those operators are applied for computing the
%     elementary operation at each step; for instance, |'roco'| means that
%     at each step, an opening by reconstruction followed by a closing by
%     reconstruction itself followed by another opening by reconstruction are 
%     applied on the input image; default: |op='roc'|.
%
% *|se|* : (optional) string or cell of defining the structuring elements used
%     for computing the granulometry; it can be either a cell of |N| structures 
%     of type |STREL|, each |se{i}| representing the SE used by the operator
%     at step |i<=N|, or a string defining the shape of the elementary SE
%     (see function |FLATSTREL|); in this latter case, the size of the SE
%     needs to be passed in the variable |s1| (see below); default: |se='disk'|.
% 
%% Property [propertyname  propertyvalues]
% *|'s1', 's2'|* : optional argument(s) further defining the series of SE's 
%     used by the elementary morphological operation at each step in the case
%     the variable |se| (see above) was passed as a string (and not a cell of
%     |STREL| structures); |s1| is a vector of size |N|, each |s1(i)| providing
%     the size of the SE used by the operator at step |i<=N|; |s2| further
%     completes the definition of the SE, depending on its type (see function
%     |FLATSTREL|); default: |s1=2:2:10| and |s2=[]|.
%   
%% Outputs
% *|G|* : a cell of size |N| representing a granulometry.
%
%% References
% [Mathe75]  G. Matheron: "Random Sets and Integral Geometry", Wiley, New
%      York (USA), 1975.
%
% [Serra82]  J. Serra: "Image Analysis and Mathematical Morphology", Academic
%      Press, New York (USA), 1982.
%
% [CD94]  Y. Chen and E.R. Dougherty: "Gray-scale morphological granulometric
%      texture classification", Opt. Eng. 33(8):2713-2722, 1994.
%      <http://spiedigitallibrary.org/oe/resource/1/opegar/v33/i8/p2713_s1>
%
% [DC99]  E.R. Dougherty and Y. Chen: "Granulometric filters", in: "Nonlinear
%      Image Filtering", E.R. Dougherty and J. Astola eds., SPIE and IEEE 
%      Presses, Belingham, pp. 121-162, 1999.
%
%% See also
% Related:
% <MORPHPROFILE.html |MORPHPROFILE|>,
% <matlab:web(whichpath('WATERSHED')) |WATERSHED|>.
% Called:
% <GRANULOMETRY_BASE.html |GRANULOMETRY_BASE|>,
% <FLATSTREL.html |FLATSTREL|>.

%% Function implementation
function [G, varargout] = granulometry( I, varargin )

%%
% parsing parameters

if isempty(ver('images'))
    error('morphprofile:errortoolbox', 'Image Processing toolbox required');
end

error(nargchk(1, 16, nargin, 'struct'));
error(nargoutchk(1, 4, nargout, 'struct'));

if ~isnumeric(I)
    error('morphprofile:inputparameter','a matrix is required in input'); 
end

p = createParser('MORPHPROFILE');   
p.addOptional('op', 'roc', @(x)ischar(x) && ...
    any(strcmpi(x,{'ro', 'ropen', 'rc', 'rclose', 'roc','rocmax', ...
    'o', 'open', 'c', 'close', 'oc', 'ocmax', ...
    'e', 'erode', 'd', 'dilate'})));
p.addOptional('se', 'disk', @(x) (iscell(x) && strcmp(class(x{1}),'strel')) || ...
    (ischar(x) && any(strcmpi(x,{'disk','rectangle','square','diamond', ...
                'line','periodicline','arbitrary','octagon','pair'}))));
p.addParamValue('s1', 2:2:10, @(x)isnumeric(x));
p.addParamValue('s2', [], @(x)isnumeric(x));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% checking/setting variables

C = size(I,3);

if ischar(p.se)
    shape = p.se;
    N = length(p.s1);
    
    if ~iscell(p.s1),
        if any(strcmp(shape,{'disk'}))
            p.s1 = floor(p.s1/2);
        end
        p.s1 = num2cell(p.s1);  
    end
    if ~iscell(p.s2),
        if isempty(p.s2)
            p.s2 = cell(N,1);  
        else
            if isscalar(p.s2),  p.s2 = repmat(p.s2, 1, N);  end;
            p.s2 = num2cell(p.s2); 
        end
    end
    
    p.se = cell(N,1);
    for i=1:N
        p.se{i} = flatstrel(shape, p.s1{i}, p.s2{i});
    end
    
end

%% 
% main processing

G = granulometry_base( I, p.op, p.se );
    
if nargout==2,  varargout{1} = p.se;  end

%%
% display

if p.disp
    
    if any(strcmp(p.op,'o')),       p.op = 'open';
    elseif any(strcmp(p.op,'c')),   p.op = 'close';
    elseif strcmp(p.op,'oc'),       p.op = 'open-close';
    elseif strcmp(p.op,'co'),       p.op = 'close-open';
    elseif strcmp(p.op,'oco'),      p.op = 'open-close-open';
    elseif strcmp(p.op,'coc'),      p.op = 'close-open-close';
    elseif any(strcmp(p.op,'ro')),  p.op = 'open-by-rec';
    elseif any(strcmp(p.op,'rc')),  p.op = 'close-by-rec';
    elseif strcmp(p.op,'roc'),      p.op = 'open-close-by-rec';
    elseif strcmp(p.op,'rco'),      p.op = 'close-open-by-rec';
    elseif strcmp(p.op,'roco'),     p.op = 'open-close-open-by-rec';
    elseif strcmp(p.op,'rcoc'),     p.op = 'close-open-close-by-rec';
    end
    
    N = length(G);
    figure; ndisp = 3; mdisp = ceil(N/ndisp);
    for i=1:N
        subplot(mdisp, ndisp, i),  imagesc(rescale(G{i})), axis image off;
    end
    if C==1,  colormap gray; end
    suptitle(['MP of ' p.op]);
end

end % end of granulometry
