%% MERGESTRUCT - Merge structures with unique fields.
%
%% Syntax
%     sout = MERGESTRUCT(varargin);
%
%% See also  
% Related:
% <CATSTRUCT.html |CATSTRUCT|>,
% <STRUCTA2ASTRUCT.html |STRUCTA2ASTRUCT|>.
% Called: 
% <matlab:webpub(whichpath('STRUCT2CELL')) |STRUCT2CELL|>,
% <matlab:webpub(whichpath('CELL2STRUCT')) |CELL2STRUCT|>,
% <matlab:webpub(whichpath('UNIQUE')) |UNIQUE|>.

%% Function implementation
function sout = mergestruct(varargin)

%%
% start with collecting fieldnames, checking implicitly
% that inputs are structures.
fn = [];
for k = 1:nargin
    try
        fn = [fn ; fieldnames(varargin{k})];                           %#ok 
    catch MEstruct
        throw(MEstruct)
    end
end

%%
% make sure the field names are unique.
if length(fn) ~= length(unique(fn))
    error('mergestruct:FieldsNotUnique',...
        'Field names must be unique');
end

%%
% now concatenate the data from each struct.  Can't use
% structfun since input structs may not be scalar.
c = [];
for k = 1:nargin
    try
        c = [c , struct2cell(varargin{k})];                            %#ok 
    catch MEdata
        throw(MEdata);
    end
end

%%
% construct the output.
sout = cell2struct(c, fn, 1);
end % end of mergestruct