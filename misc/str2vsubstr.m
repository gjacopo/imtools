%% STR2VSUBSTR - Concatenate strings vertically.
%
%% Description
% Vertically concatenate the substrings of a collection of patterns strings
% found in the input string.
% 
%% Syntax
%      [vstr, nsub] = STR2VSUBSTR(str, subvstr);
%
%% Inputs
% *|str|* : a string with concatenated substrings to extract.
%
% *|subvstr|* : an array of all the possible substring patterns to look for;
%     the substring patterns are stored in the lines of this array (see
%     function |STRVCAT|).
%
%% Outputs
% *|vstr|* : the substring patterns of subvstr found in |str| rewritten by
%     vertically concatenating those substrings in the order of their 
%     appearance in |str|.
%
% *|nsub|* : optional output storing the number of substrings found.
%
%% See also  
% Called: 
% <matlab:webpub(whichpath('STRTRIM')) |STRTRIM|>,
% <matlab:webpub(whichpath('STRFIND')) |STRFIND|>,
% <matlab:webpub(whichpath('SORT')) |SORT|>.

%% Function implementation
function [vstr, nsub] = str2vsubstr(str, subvstr)

error(nargchk(1, 2, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

vstr = subvstr;
nsub = length(vstr); % number of computed features
isub = []; % stores the order of appearance of the features in the
% string feature, which will also be the order of the feature in the
% output list of arguments

i=1;
while i<=nsub
    pattern = strtrim(vstr(i,:));
    j = strfind(lower(str),pattern);
    if isempty(j),   vstr(i,:) = []; % the feature wont be estimated
    else
        isub = [isub j];    %#ok   position of appearance of the pattern 
        str(j+(0:length(pattern)-1)) = ''; % delete the pattern from str
        i = i+1; 
    end;
    nsub = size(vstr,1);
end

[~,j] = sort(isub); % order the appearance in the string
vstr = vstr(j,:); % reorder the features for the output
end % end of str2vsubstr