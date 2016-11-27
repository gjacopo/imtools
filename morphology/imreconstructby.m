%% IMRECONSTRUCTBY - Opening- and Closing-by-Reconstruction.
%
%% Description
% Perform opening- and closing-by-reconstruction, using the function
% |IMRECONSTRUCT| of the Image Processing toolbox.
%
%% Syntax
%     R = IMRECONSTRUCTBY(I);
%     R = IMRECONSTRUCTBY(I, op);
%     R = IMRECONSTRUCTBY(I, op, shape);
%     R = IMRECONSTRUCTBY(I, op, shape, s1);
%     R = IMRECONSTRUCTBY(I, op, shape, s1, s2);
%
%% Inputs
% *|I|* : imput image.
%
% *|op|* : string specifying if the opening (|'ro'| or |'ropen'|) or the
%     closing (|'rc'| or |'rclose'|) by reconstruction is applied; default:
%     |op='ro'|.
%
% *|shape, s1, s2|* : optional arguments defining the SE used by the operator;
%     see function |FLATSTREL|; default: |shape='disk', s1=3:2:10, s2=[]|.
%
%% Outputs
% *|R|* : filtered image.
%
%% See also
% Related:
% <matlab:web(whichpath('IMRECONSTRUCT')) |IMRECONSTRUCT|>,
% <IMROPEN.html |IMROPEN|>,
% <IMRCLOSE.html |IMRCLOSE|>,
% <IMRECONSTRUCTBY.html |IMRECONSTRUCTBY|>.
% Called: 
% <IMRECONSTRUCTBY_BASE.html |IMRECONSTRUCTBY_BASE|>.
% <FLATSTREL.html |FLATSTREL|>.

%% Function implementation
function R = imreconstructby(I, varargin)

if isempty(ver('images'))
    error('imreconstructby:errortoolbox', 'Image Processing toolbox required');
end

error(nargchk(1, 13, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

%% Parsing parameters

if ~isnumeric(I)
    error('imreconstructby:inputparameter','a matrix is required in input'); 
end

% optional parameters
p = createParser('IMRECONSTRUCTBY');   

p.addOptional('op', 'ro', @(x)ischar(x) && any(strcmpi(x,{'ro','ropen','rc','rclose'})));
p.addOptional('se', 'disk', @(x) strcmp(class(x),'strel') || ...
    (ischar(x) && any(strcmpi(x,{'disk','rectangle','square','diamond', ...
                'line','periodicline','arbitrary','octagon','pair'}))));
p.addOptional('s1', 3:2:10, @(x)isnumeric(x));
p.addOptional('s2', [], @(x)isnumeric(x));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            


%% Checking/setting variables

if ischar(p.se)
    shape = p.se;
    N = length(p.s1);
    
    if ~iscell(p.s1),  p.s1 = num2cell(p.s1);  end
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

R = imreconstructby_base(I, p.op, p.se);

%%
% display

if p.disp
    if any(strcmp(p.op,{'ro','ropen'})),       p.op = 'open-by-rec';
    elseif any(strcmp(p.op,{'rc','rclose'})),  p.op = 'close-by-rec';
    end

    figure, imagesc(rescale(R)), axis image off, title(p.op);
    if size(R,3)==1,  colormap gray;  end
end

end % end of imreconstructby
