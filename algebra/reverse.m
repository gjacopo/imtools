%% REVERSE - Flip a vector.
%
%% Syntax
%   r = REVERSE(x);
%
%% Remark
% Useless, use |FLIPLR| instead !!!
%
%% See also  
% Related:
% <RESCALE.html RESCALE>.

%% Function implementation
function r = reverse(x)
r = x(end:-1:1);
end