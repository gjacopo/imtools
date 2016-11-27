%% FINDLOCALMAX - Local maxima of a matrix.
%
%% Description
% Determine the local maxima of a given matrix image.
%
%% Syntax
%     max_local = FINDLOCALMAX(I);  
%     [max_local,row,col] = FINDLOCALMAX(I, radius, method );  
%
%% Inputs
% *|I|* : a |(N,M)| matrix containing values to filter.
%
% *|radius : optional variable setting the radius of the neighborhood; 
%     default: |radius=1|, ie maxima in neighbourhoods |(3,3)| are extracted.
%
% *|method|* : optional string defining the technique used for extracting
%     the maxima; it is either |'dil'| (fast, not unique), |'filt'| (fast) or
%     'screen' (slow); default: |method='filt'|.
%
%% Outputs
% *|row|* : the row position of the local maxima.
%
% *|col|* : the column position of the local maxima.
%
% *|max_local|*	: the |(N,M)| matrix containing values of |val| on unique
%       local maximum.
%
%% Remark
% |FINDLOCALMAX| and |FINDLOCALEXTREMA(..,'ord','max','corr','forall')|
% provide identical results, except for the borders of the image.
%
%% See also
% Related:
% <FINDLOCALEXTREMA.html |FINDLOCALEXTREMA|>,
% <../../sharpen/html/MAPTRANSITION.html |MAPTRANSITION|>,
% <matlab:webpub(whichpath('ORDFILT2')) |ORDFILT2|>.
% Called: 
% <FINDLOCALMAX_BASE.html |FINDLOCALMAX_BASE|>.

%% Function implementation
function [max_local,row,col] = findlocalmax(I,varargin)

%%
% parsing and checking parameters
if ~isnumeric(I)
    error('findlocalmax:inputerror','a matrix is required in input'); 
end

p = createParser('FINDLOCALMAX');   
% mandatory parameter
% optional parameters
p.addOptional('radius',1., @(x)x>=0); % just for testing
p.addOptional('find', 'filt', @(x)ischar(x) && ...
    any(strcmpi(x,{'filt','dil','screen'})));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%% 
% dealing with multispectral images
% [M,N,C] = size(I);
% if C>1
%     max_local = zeros(size(I));
%     row = cell(C); col = cell(C);
%     for c=1:C
%         [max_local(:,:,c), row{c}, col{c}] = findlocalmax(I(:,:,c), 'parse', p);
%     end
%     return;
% end

%% 
% main computation

[max_local,row,col] = findlocalmax_base(I, p.radius, p.find);

%%
% display

if p.disp
    figure, imagesc(max_local), axis image off, title('local maxima')
end

end % end of findlocalmax
