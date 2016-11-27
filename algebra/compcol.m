%% COMPCOL - Compare columns of a matrix.
%
%% Description
% Compare the columns of |A| in increasing order (ie, |A(:,1)| 'is compared 
% to' |A(:,2), A(:,2)| 'is compared to' |A(:,3),...|) where the comparison
% function 'is compared to' is given by the function handled by the string
% |order|.
%
%% Syntax
%       R = COMPCOL(A);
%       R = COMPCOL(A,order);
%
%% Inputs
% *|A|* : input matrix whose colums have to be compared; if |dims(A)==3|,
%     then the matrix |A| is reshaped so that its elements are compared along
%     the 3rd dimension.
%
% *|order|* : optional input string handling the function used for comparing
%     the columns of |A|; order is either |'>'|, |'<'|, |'<='|, |'>='| or 
%     |'=='|; default: vorder='<'|.
% 
%% Output
% *|R|* : a boolean matrix with same number of lines as |A| stating for each
%     line |i| of the matrix |A| if the assertion:
%            |A(i,1) order A(i,2) order ... order A(i,ncols)|
%     is |true| or |false|.
%
%% See also
% Related: 
% <matlab:webpub(whichpath('SORT')) |SORT|>,
% <matlab:webpub(whichpath('SORTROWS')) |SORTROWS|>.

%% Function implementation
function R = compcol(A,varargin)

%%
% parsing the parameters

if ~isnumeric(A)
    error('compcol:inputerror',...
        'a matrix is required in input'); 
end

p = createParser('COMPCOL');   % create an instance of the inputParser class.
p.addOptional('order', '<', @(x)ischar(x) && ...
    any(strcmpi(x,{'>','<','<=','>=','=='})));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            


%%
% internal variables

ndims = nb_dims(A);
[x y ncols] = size(A);
if ndims==2 || ncols==1
    ncols = y;
elseif ndims==3
    A = reshape(A,[x*y ncols]);
    x = x*y;
elseif ndims>3
    error('compcol:inputerror',...
        'cannot process vectors or matrices with dimension>3');
end

% set the order variable to the handle for the comparison function selected
order = str2func(p.order);

%% 
% testing

% initialize the output image
R = true(x,1); % ones(x,1);

%%
% main comparison function
for i=1:ncols-1
    % R = eval(['R & A(:,i)' order 'A(:,i+1)']);
    R = R & order(A(:,i),A(:,i+1));
    % example: if order='<' and A with size ncols=3, this expresion is then 
    % equivalent to:
    %     R = A(:,1)<I(:,2) & A(:,2)<I(:3)
end

if ndims==3
    A = reshape(A,[x/y y ncols]);                                      %#ok
    R = reshape(R,[x/y y]);
end

end% end of compcol
