%% FLATSTREL - Wrapper utility for function |STREL|.
%
%% Description
% Nothing else than a wrapper utility function for |STREL| used for creating
% 2D flat structuring elements (SE). All types of SE's are supported,
% except for |'ball'|.
%
%% Syntax
%     S = FLATSTREL(shape);
%     S = FLATSTREL(s1);
%     S = FLATSTREL(shape, s1);
%     S = FLATSTREL(s1, shape);
%     S = FLATSTREL(shape, s1, s2);
%
%% Inputs
% *|shape|* : any string among |'disk'|, |'rectangle'|, |'square'|, |'diamond'|,
%     |'line'|, |'periodicline'|, |'arbitrary'|, |'octagon'| and |'pair'|
%     (note that |'ball'| is not supported) or any scalar value. In the former
%     case, a default SE of type |shape| and size 3 will be output, in the latter
%     case, a SE of type |'disk'| and size shape will be output.
%
% *|s1|* : optional shape (string, see above) or size of the SE when |nargin=2|;
%     when |nargin=3| (ie. |s2| is also defined), |s1| is necessarly a scalar
%     value.
%
% *|s2|* : optional scalar value for defining the SE (see |STREL| possibly defined
%     with 2 parameters: |'periodicline'|, |'disk'|, |'line'|).
%
%% Outputs
% *|S|* : structuring element of class |STREL|.
%
%% See also
% Called: 
% <matlab:web(whichpath('STREL')) |STREL|>.
%% Function implementation
function S = flatstrel(shape, s1, s2)

if isempty(ver('images'))
    error('flatstrel:errortoolbox', 'Image Processing toolbox required');
end

error(nargchk(1, 3, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

if nargin==1 || ((nargin==2 && isempty(s1)) || (nargin==3 && isempty(s1) && isempty(s2)))
    if isnumeric(shape),  S = strel('disk', shape, 0);
    elseif ischar(shape)
        if ~any(strcmp(shape,{'disk','rectangle','square','diamond', ...
                'line','periodicline','arbitrary','octagon','pair'}))
            error('flatstrel:errorinput', ['shape ''' shape ''' of the SE not supported']);
        end
        if strcmp(shape,'line'),      S = strel(shape, 3, 0);
        elseif strcmp(shape,'disk'),  S = strel(shape, 3, 0);
        else                          S = strel(shape, 3);
        end
    else
        error('flatstrel:errorinput', ...
            'default size or shape of the SE required in input');
    end
    
elseif nargin==2 || (nargin==3 && isempty(s2))
    if isnumeric(shape) && ischar(s1)
        tmp = shape; shape = s1; s1 = tmp;
    end
    if ischar(shape)
        if ~any(strcmp(shape,{'disk','rectangle','square','diamond', ...
                'line','periodicline','arbitrary','octagon','pair'}))
            error('flatstrel:errorinput', ['shape ''' shape ''' of the SE not supported']);
        end
        if  isnumeric(s1)
            if strcmp(shape,'line'),              S = strel(shape, s1, 0);
            elseif strcmp(shape,'disk'),          S = strel(shape, s1, 4);
            elseif strcmp(shape,'periodicline'),  S = strel(shape, s1, 2);
            elseif strcmp(shape,'rectangle'),     S = strel(shape, [s1, s1]);
            else                                  S = strel(shape, s1);
            end
        end
    else
        error('flatstrel:errorinput', 'size and shape of the SE required in input');
    end
        
else
    if ~(ischar(shape) && isnumeric(s1) && isnumeric(s2))
        error('flatstrel:errorinput', ...
            'shape and size parameters of the SE required in this order in input');
    elseif ~any(strcmp(shape,{'disk','rectangle','square','diamond', ...
            'line','periodicline','arbitrary','octagon','pair'}))
        error('flatstrel:errorinput', ['shape ''' shape ''' of the SE not supported']);
    elseif any(strcmp(shape,{'disk','line','arbitrary'})),  S = strel(shape, s1, s2);
    elseif strcmp(shape,'rectangle'),                       S = strel(shape, [s1, s2]);
    else                                                    S = strel(shape, s1);
    end
end

end % end of flatstrel
