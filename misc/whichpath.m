%% WHICHPATH - Function full path name.
%
%% Description
% Similarly to WHICH (which it uses), return the full path name of a function,
% but filter out of the extra information output by that former function.
%
%% Syntax
%       p = WHICHPATH(filename);   
%
%% See also
% Called:
% <matlab:webpub(whichpath('WHICH')) |WHICH|>,
% <matlab:webpub(whichpath('STRFIND')) |STRFIND|>.

%% Function implementation
function p = whichpath(filename)

if ~(exist(filename,'file') || exist(filename,'builtin'))
    p = [];
    return
end

if ~strcmpi(filename(end-1:end),'.m'),  filename = [filename '.m'];  end

% [~, p] = unix(['find ' matlabroot ' -name ' filename]);

p = which(filename);
i = strfind(p,'(');  j = strfind(p,')');
if ~isempty(i) && ~isempty(j)
  p = p(i+1:j-1);
end

end % end of whichpath
