%% NB_DIMS - Overload NDIMS ('debugged' version).
%
%% Syntax
%     d = NB_DIMS(x);
%
%% See also
% Related:
% <matlab:webpub(whichpath('NDIMS')) |NDIMS|>.

%% Function implementation
function d = nb_dims(x)

if isempty(x)
    d = 0;
    return;
end

d = ndims(x);

if d==2 && (size(x,1)==1 || size(x,2)==1)
    d = 1;
end % end of nb_dims