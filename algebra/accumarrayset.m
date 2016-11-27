%% ACCUMARRAYSET - Overload |ACCUMARRAY|.
%
%% Description
% Group elements from a data vector set as they share a same index. Fast
% version of |ACCUMARRAY(subs, val, [], @(x){x})|.
%
%% Syntax
%     A = accumarrayset(subs);
%     [A, T] = accumarrayset(subs, val);
%     [A, T] = accumarrayset(subs, val, spar);
%
%% Inputs
% *|subs, val|* : the position of an element in the vector |subs| determines
%     which value of the vector |vals| it selects for the accumulated vector;
%     the value of an element in |subs| determines the position of the
%     accumulated vector in the output; see function |ACCUMARRAY|; |val| is
%     optional, in the case it is not passed the vector |(1:numel(subs))| is
%     used by default.
%
% *|spar|* : (optional) logical flag defining if the output cell array |A| 
%     (see below) should contain empty cells; when set to |true|, only those 
%     indices found in |subs| are represented in the final output; which ones
%     are represented is given in the vector |T|; default: |spar=true|.
%
%% Outputs
% *|A|* : cell array grouping values of |val| with same index in |subs|.
%
% *|T|* : sparse matrix giving for each (unique) value/entry found in |subs|
%     the position (index) in the cell |A| of the set grouping the corresponding 
%     elements in |val|; |T| is output when the variable |spar| is set to 
%     |true|. 
%
%% See also
% Related:
% <matlab:webpub(whichpath('ACCUMARRAY')) |ACCUMARRAY|>.
% Called:
% <MAT2RC.html |MAT2RC|>.
% <matlab:webpub(whichpath('SORT')) |SORT|>,
% <matlab:webpub(whichpath('FIND')) |FIND|>,
% <matlab:webpub(whichpath('DIFF')) |DIFF|>,
% <matlab:webpub(whichpath('SPARSE')) |SPARSE|>.

%% Function implementation
function [A, T] = accumarrayset(subs, val, spar)

if nargin<3,  spar = false;
    if nargin<2,  val = [];  end
end

if isempty(val),  val = (1:numel(subs))';  end
subs = mat2rc(subs,'c');  val = mat2rc(val,'c');
% tic
% A = accumarray(subs, val, [], @(x){x});
% toc

%tic
[r,c] = size(subs);
L = zeros(r,1);
L(:) = val; % (1:r);
L = L(:,ones(1,c));
[T,idx] = sort(subs(:)); % we need T for finding rows
L = L(idx);
idx = [0; find(diff(T)); numel(T)]; % run length encoding
T = T(idx(2:end));
t = numel(T);
if spar
    A = cell(t, 1);
    for ii = 1:t
        A{ii} = L(idx(ii)+1:idx(ii+1));
    end
    T = sparse(T,ones(size(T)),1:t,T(end),1);
else
    A = cell(T(end), 1);
    % create same look up table as ACCUMARRAY does
    for ii = 1:t
        A{T(ii)} = L(idx(ii)+1:idx(ii+1));
    end
end
% toc

end % end of accumarrayset