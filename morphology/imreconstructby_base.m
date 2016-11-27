%% IMRECONSTRUCTBY_BASE - Base function for IMRECONSTRUCTBY.
%
%% Syntax
%     R = IMRECONSTRUCTBY_BASE(I, op);
%     R = IMRECONSTRUCTBY_BASE(I, op, shape);
%     R = IMRECONSTRUCTBY_BASE(I, op, shape, s1);
%     R = IMRECONSTRUCTBY_BASE(I, op, shape, s1, s2);
%
%% See also
% Related:
% <matlab:web(whichpath('IMRECONSTRUCT')) |IMRECONSTRUCT|>,
% <IMRECONSTRUCTBY.html |IMRECONSTRUCTBY|>.
% Called: 
% <matlab:web(whichpath('IMRECONSTRUCT')) |IMRECONSTRUCT|>,
% <matlab:web(whichpath('IMERODE')) |IMERODE|>,
% <matlab:web(whichpath('IMDILATE')) |IMDILATE|>,
% <FLATSTREL.html |FLATSTREL|>.

%% Function implementation
function R = imreconstructby_base(I, op, varargin)

% error(nargchk(1, 5, nargin, 'struct'));
% error(nargoutchk(1, 1, nargout, 'struct'));

%%
% preparing parameters

if nargin<=2,  
    se = strel('disk',1);
    if nargin==1,  op = 'ropen'; end

elseif nargin==3 && strcmp(class(varargin{1}),'strel')
    se = varargin{1};

else
    se = flatstrel(varargin{:});
end

% %% Dealing with multispectral input !!! NO NEED!!!
% [X,Y,C] = size(I); % possibly multichannel when C>1 
% if C>1
%     R = zeros(size(I));
%     for c=1:C
%         R(:,:,c) = imreconstructby_base(I(:,:,c), op, se);
%     end
% end

%% 
% main computation

switch op
    case {'ro','ropen'} 
        % opening-by-reconstruction using imerode followed by imreconstruct
        Ie = imerode(I, se); % marker
        R = imreconstruct(Ie, I);
        
    case {'rc','rclose'}
        % closing-by-reconstruction using imdilate followed by imreconstruct
        % inputs and output of imreconstruct must be complemented
        Id = imcomplement(imdilate(I, se)); % marker
        R = imcomplement(imreconstruct(Id, imcomplement(I)));
        % it is also:
        % R2 = imcomplement(imreconstructby_base(imcomplement(I), 'ro', se));
        
end

end % end of imreconstructby_base
