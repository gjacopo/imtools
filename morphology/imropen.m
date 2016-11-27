%% IMROPEN - Opening-by-Reconstruction.
%
%% Description
% Wrapping function for performing opening by reconstruction (see function
% |IMRECONSTRUCTBY_BASE|).
%
%% Syntax
%     F = IMROPEN(I, se);
%     F = IMROPEN(I, nhood);
%
%% Inputs
% *|I|* : imput image.
%
% *|se|* : opening will be performed with the single structuring element (SE)
%     defined in |se|.
%
% *|nhood|* : opening will be performed with the SE |strel(nhood)|, where
%     |nhood| is an array of 0's and 1's that specifies the SE neighborhood.
%
%% Outputs
% *|F|* : filtered image.
%
%% See also
% Related:
% <matlab:web(whichpath('IMOPEN')) |IMOPEN|>,
% <matlab:web(whichpath('IMRECONSTRUCT')) |IMRECONSTRUCT|>,
% <IMRCLOSE.html |IMRCLOSE|>,
% <IMRECONSTRUCTBY.html |IMRECONSTRUCTBY|>.
% Called: 
% <IMRECONSTRUCTBY_BASE.html |IMRECONSTRUCTBY_BASE|>.

%% Function implementation
function F = imropen(I, varargin)

error(nargchk(1, 10, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

%%
% parsing parameters

if ~isnumeric(I)
    error('imropen:inputparameter','a matrix is required in input'); 
end

p = createParser('IMROPEN');   
p.addRequired('se', @(x)strcmp(class(x),'strel') || ...
    (isnumeric(x) && (islogical(x) || all(ismember(x(:),[0,1])))));
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% checking/setting variables

if ~isequal(class(p.se),'strel')
    p.se = strel(p.se);
end

F = imreconstructby_base(I, 'ropen', p.se);

if p.disp
   figure, imagesc(rescale(F)), axis image off, title('opening by reconstruction');
   if size(I,3)==1,  colormap gray;  end
end

end % end of imropen
