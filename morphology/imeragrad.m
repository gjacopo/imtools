%% IMERAGRAD - Morphological gradient ERAGRAD.
%
%% Description
% Compute the so-called ERAGRAD operation described in [BG09] consisting of
% a morphological gradient followed by an erosion.
%
%% Syntax
%      E = IMERAGRAD(I);
%      E = IMERAGRAD(I, shape);
%      E = IMERAGRAD(I, shape, s1);
%      E = IMERAGRAD(I, shape, s1, s2);
%
%% Inputs
% *|I|* : input image of size |(X,Y,C)|, with |C>1| when |I| is multispectral.
%
% *|shape, s1, s2|* : optional arguments defining the SE used by the elementary
%     operators called by |IMERAGRAD|; see function |FLATSTREL| for further
%     explanation; default: |shape='square', s1=3, s2=[]|.
%
%% Outputs
% *|E|* : filtered image of same size as |I| output by the ERAGRAD operator.
%
%% Reference
% [BG09]  R. Bellens and S. Gautama: "Erosion after gradient (ErAGrad)
%      morphological profile", Proc. IEEE International Geoscience & Remote
%      Sensing Symposium, 2009.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5418079>
%
%% See also
% Related:
% <matlab:web(whichpath('IMCLOSE')) |IMCLOSE|>,
% <IMRCLOSE.html |IMRCLOSE|>,
% <IMRECONSTRUCTBY.html |IMRECONSTRUCTBY|>.
% Called:
% <IMERAGRAD_BASE.html |IMERAGRAD_BASE|>,
% <FLATSTREL.html |FLATSTREL|>.

%% Function implementation
function E = imeragrad(I, varargin)

if isempty(ver('images'))
    error('imeragrad:errortoolbox', 'Image Processing toolbox required');
end

error(nargchk(1, 12, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

%%
% parsing parameters

if ~isnumeric(I)
    error('imeragrad:inputparameter','a matrix is required in input'); 
end

% optional parameters
p = createParser('IMERAGRAD');   

p.addOptional('se', 'disk', @(x) strcmp(class(x),'strel') || ...
    (ischar(x) && any(strcmpi(x,{'disk','rectangle','square','diamond', ...
                'line','periodicline','arbitrary','octagon','pair'}))));
p.addOptional('s1', 3:2:10, @(x)isnumeric(x));
p.addOptional('s2', [], @(x)isnumeric(x));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
%  checking/setting variables

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

E = imeragrad_base(I,p.se);

%%
% display

if p.disp
    figure, imagesc(rescale(E)), axis image off, title('erosion after gradient');
    if size(E,3)==1,  colormap gray;  end
end

end % end of imeragrad
