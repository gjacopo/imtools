%% IMERAGRAD_BASE - Base function for IMERAGRAD.
%
%% Syntax
%     E = IMERAGRAD_BASE(I);
%     E = IMERAGRAD_BASE(I, shape);
%     E = IMERAGRAD_BASE(I, shape, s1);
%     E = IMERAGRAD_BASE(I, shape, s1, s2);
%
%% See also
% Related:
% <IMERAGRAD.html |IMERAGRAD|>,
% <IMRECONSTRUCTBY_BASE.html |IMRECONSTRUCTBY_BASE|>.
% Called:
% <matlab:web(whichpath('IMDILATE')) |IMDILATE|>,
% <matlab:web(whichpath('IMERODE')) |IMERODE|>.

%% Function implementation
function E = imeragrad_base(I, varargin)

% error(nargchk(1, 4, nargin, 'struct'));
% error(nargoutchk(1, 1, nargout, 'struct'));

if nargin<=2,  
    se = strel('square',3);
elseif nargin==3 && strcmp(class(varargin{1}),'strel')
    se = varargin{1};
else
    se = flatstrel(varargin{:});
end

E = imdilate(I, se) - imerode(I, se); % morphological gradient
E = imerode(E, se);

end % end of imeragrad_base