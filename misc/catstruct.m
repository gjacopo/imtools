%% CATSTRUCT - Concatenate structures.
%
%% Syntax _(i)_
%     X = CATSTRUCT(S1,S2,S3,...) 
%
% Concate the structures |S1|, |S2|, ... into one structure |X|. 
%
%% Syntax _(ii)_ 
%     CATSTRUCT(S1,S2,'sorted') 
%
% Sort the fieldnames alphabetically.
%
%% Remarks
% * If a fieldname occurs more than once in the argument list, only the last
%   occurence is used, and the fields are alphabetically sorted.
%
% * To sort the fieldnames of a structure |A| use:
%
%     A = CATSTRUCT(A,'sorted') ;
%
% * To concatenate two similar array of structs use simple concatenation:
%
%     A = dir('*.mat') ; B = dir('*.m') ; C = [A ; B] ;
%
% * When there is nothing to concatenate, the result will be an empty
%   struct (|(0,0)| struct array with no fields). 
%
%% Example:
%     A.name = 'Me' ; B.income = 99999 ; 
%     X = CATSTRUCT(A,B) 
%     % -> X.name = 'Me';  X.income = 99999 ;
%
%% Acknowledgment
% for Matlab R13 and up - version 2.2 (oct 2008)
% <mailto:jos@jasen.nl Jos van der Geest>
%
%% See also
% Related:
% <MERGESTRUCT.html |MERGESTRUCT|>,
% <STRUCTA2ASTRUCT.html |STRUCTA2ASTRUCT|>.
% Called:
% <matlab:webpub(whichpath('CAT')) |CAT|>,
% <matlab:webpub(whichpath('STRUCT')) |STRUCT|>,
% <matlab:webpub(whichpath('FIELDNAMES')) |FIELDNAMES|>,
% <matlab:webpub(whichpath('STRUCT2CELL')) |STRUCT2CELL|>.

% History
% Created:  2005
% Revisions
%   2.0 (sep 2007) removed bug when dealing with fields containing cell
%                  arrays (Thanks to Rene Willemink) 
%   2.1 (sep 2008) added warning and error identifiers
%   2.2 (oct 2008) fixed error when dealing with empty structs (Thanks to
%                  Lars Barring)

%% Function implementation
function A = catstruct(varargin)

N = nargin ;

error(nargchk(1,Inf,N)) ;

if ~isstruct(varargin{end}),
    if isequal(varargin{end},'sorted'),
        sorted = 1 ;
        N = N-1 ;
        if N < 1,
            A = struct([]) ;
            return
        end
    else
        error('catstruct:InvalidArgument','Last argument should be a structure, or the string "sorted".') ;
    end
else
    sorted = 0 ;
end

FN = cell(N,1) ;
VAL = cell(N,1) ;

for ii=1:N,
    X = varargin{ii} ;
    if ~isstruct(X),
        error('catstruct:InvalidArgument',['Argument #' num2str(ii) ' is not a structure.']) ;
    end
    if ~isempty(X),
        % empty structs are ignored
        FN{ii} = fieldnames(X) ;
        VAL{ii} = struct2cell(X) ; 
    end
end

FN = cat(1,FN{:}) ;
VAL = cat(1,VAL{:}) ;
[UFN,ind] = unique(FN) ;

if numel(UFN) ~= numel(FN),
    warning('catstruct:DuplicatesFound','Duplicate fieldnames found. Last value is used and fields are sorted') ;
    sorted = 1 ;
end

if sorted,
    VAL = VAL(ind) ;
    FN = FN(ind) ;
end

if ~isempty(FN),
    % This deals correctly with cell arrays
    A = cell2struct(VAL, FN);
else
    A = struct([]) ;
end % end of catstruct




