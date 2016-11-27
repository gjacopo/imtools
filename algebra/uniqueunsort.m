%% UNIQUEUNSORT - Unsorted list of unique elements.
%
%% Description
% Find unique elements of vector likewise |UNIQUE|, but, contrary to it,
% without sorting the results.
%
%% Syntax
%       Y = UNIQUEUNSORT(X);
%       [Y,U] = UNIQUEUNSORT(X);
%       [Y,U] = UNIQUEUNSORT(X,'first'/'last');
%       [Y,U] = UNIQUEUNSORT(X,'rows');
%       [Y,U] = UNIQUEUNSORT(X,'first'/'last','rows');
% 
%% Inputs
% *|X|* : vector or matrix where to extract unique values.
%
% *|'first', 'last', 'rows'|* : see |UNIQUE| help; default: |'first'| only
%      is set (when using |UNIQUE|, |'last'| is the default setting).
% 
%% Output
% *|Y|* : unique values contained in the input data, without sorting, i.e. 
%     output as the way they are 'read' in the input |X| - from the beginning
%     in the case |'first'|, from the end in the case |'last'|.
%
%% See also
% Related:
% <matlab:webpub(whichpath('UNIQUE')) |UNIQUE|>.
% Called:
% <matlab:webpub(whichpath('SORT')) |SORT|>,
% <matlab:webpub(whichpath('SORTROWS')) |SORTROWS|>,
% <matlab:webpub(whichpath('DIFF')) |DIFF|>.

%% Function implementation
function [Y,U] = uniqueunsort(X, varargin)

%%
% check/set input variables

rowsort = [];
order = [];

for i=1:nargin-1
    if ~(ischar(varargin{i}) && any(strcmpi(varargin{i},{'first','last','rows'})))
        error('uniqueunsort:errorinput', ['unknown method ' varargin{i}])
    else
        switch varargin{i}
            case {'first','last'}
                if isempty(order),  order = varargin{i};
                else  error('uniqueunsort:warninginput', ...
                        'repeated entries ignored');
                end
                
            case 'rows'
                if isempty(rowsort),  rowsort = true;
                else  error('uniqueunsort:warninginput', ...
                        'repeated entries ignored');
                end
        end
    end
end

if isempty(order),  order = 'first';  end
if isempty(rowsort),  rowsort = false;  end

if rowsort
    fsort = @(x) sortrows(x);
else
    fsort = @sort;
    X = X(:);
end

%%
% sort the data first, retaining the indices from the sorting operation
[Xs,sortX] = fsort(X);

%%
% take the difference of the sorted results and find where the differences
% are not zero (i.e., they are different values): create the correct indices
% for these now unique values in the logical vector U
U = false(size(X,1),1);
if strcmpi(order,'last'),  Xs = [sum(diff(Xs,1),2); 1]~=0;
else                       Xs = [1; sum(diff(Xs,1),2)]~=0;
end
U(sortX) = Xs;

%%
% use this set of logical indices to extract the required values from the
% original data
Y = X(U,:);

end % end of uniqueunsort
