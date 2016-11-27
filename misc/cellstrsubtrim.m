%% CELLSTRSUBTRIM - Trim a set of strings from its matching substrings. 
%
%% Description
% Discard all the strings of a set (array or cell) of strings that match
% substrings of longer strings in the same set.
%
%% Syntax
%       [C, I] = CELLSTRSUBTRIM(S);
%       [C, I] = CELLSTRSUBTRIM(S, traversal, despace, sensitive);
%
%% Inputs
% *|S|* : cell or array of strings (possibly containing whitespace characters
%     |' '| that can be ignored, see |despace| below) from which redundant
%     elements matching substrings of longer strings are removed.
%
% *|traversal|* : (optional) flag setting the rule for defining substrings;
%     a string in |S| is considered as a substring if it matches a substring
%     of a longer string of |S|:
%  
% * when reading both of them *from left to right* in the case |traversal=1|,          
% * when reading both of them *from right to left* in the case |traversal=-1|,  
% * *indepently of the order these strings are read* in the case |traversal=0|;     
%          
% default: |traversal=0|.    
%
% *|despace|* : (optional) logical flag stating if whitespace characters |' '|
%     should be removed from the strings (where |ISSPACE| returns |true|);
%     note moreover that leading and trailing whitespace characters are 
%     automatically removed (using |STRTRIM|); default: |despace=false|.
%
% *|sensitive|* : (optional) boolean flag defining is the comparison is case
%     sensitive (|true|, then use |STRCNCMP|) or not (|false|, then use 
%     |STRNCMPI| instead); default: |sensitive=true|.
%
%% Outputs
% *|C|* : output cell (or array, depending on the class of the input |S|)
%     of strings where the substrings of |S| have been removed following the
%     rule defined by |traversal| and |despace|; note that |C| is ordered 
%     according to the length of its strings.
%     
% *|I|* : indices (referring to the input |S|) of the strings present in the 
%     output.
%     
%% Example
%  a=['fout  ';'piut  ';'pout  ';'meu g ';'meugh ';'tuipl ';  'foutre'];
%  [c,i] = cellstrsubtrim(a)
%  [c,i] = cellstrsubtrim(a,1)
%  [c,i] = cellstrsubtrim(a,-1,true)
%  [c,i] = cellstrsubtrim(a,0,true)
%  a={ ['fout'] ['piut  ']' ['po ut'] ['meu g '] ['   meugh ']' ['tuipl '] ['fou  tre']}
%  [c,i] = cellstrsubtrim(a,1,true)
%  [c,i] = cellstrsubtrim(a,0,true)
%  % both examples should return identical outputs when setting traversal=0 
%  % and despace=true
%
%% See also
% Related:
% <../../algebra/html/CELLNUMSUBTRIM.html |CELLNUMSUBTRIM|>.
% Called:
% <matlab:webpub(whichpath('TRIL')) |TRIL|>,
% <matlab:webpub(whichpath('STRNCMP')) |STRNCMP|>,
% <matlab:webpub(whichpath('STRNCMPI')) |STRNCMPI|>,
% <matlab:webpub(whichpath('STRTRIM')) |STRTRIM|>,
% <matlab:webpub(whichpath('ISSPACE')) |ISSPACE|>,
% <matlab:webpub(whichpath('CELLSTR')) |CELLSTR|>,
% <matlab:webpub(whichpath('CELL2MAT')) |CELL2MAT|>,
% <matlab:webpub(whichpath('CELLFUN')) |CELLFUN|>,
% <matlab:webpub(whichpath('RESHAPE')) |RESHAPE|>.

%% Function implementation
%--------------------------------------------------------------------------
function [C, I] = cellstrsubtrim(C, traversal, despace, sensitive)

%%
% check/set variables

if ~iscellstr(C)
    if ~ischar(C)
        error('cellstrsubtrim:errorinput', ...
            'input variable ''C'' must be a string array or a cell of strings');
    end
    isarray = true;
    C = cellstr(C);
else 
    isarray = false;
end

if nargin<4,  sensitive = true;
    if nargin<3,  despace = false;
        if nargin<2,  traversal = 0;  end
    end
end

if ~(islogical(despace) && islogical(sensitive))
    error('cellstrsubtrim:errorinput', ...
        'input variables ''sensitive'' and ''despace'' must be logical');
elseif ~isscalar(traversal) || ~ismember(traversal,[0 1 -1])
    error('cellstrsubtrim:errorinput', ...
        'input variable ''traversal'' takes values in {-1,0,1}');
end

nC = numel(C); 
ind = (1:nC);

if nC==1 % if one string only, nothing to do
    if isarray,  C = char(C);  end
    return;
end

if sensitive,  strcompare = @strncmp;
else           strcompare = @strncmpi;
end

%%
% we need to ensure that the strings are 'row strings'...
iscol = cellfun(@(c) size(c,1)>1, C);
if any(iscol),  C(iscol) = cellfun(@transpose, C(iscol), 'Uniform', false);  end

%%
% clean up a bit: first deblank (remove leading and trailing whitespace
% characters from string, not optional), then possibly'despace' (trim all
% internal whitespaces)
C = cellfun(@strtrim, C, 'Uniform', false); 
if despace,  C = cellfun(@(c) c(~isspace(c)), C, 'Uniform', false);  end

%%
% in general, get rid of empty strings
I = cellfun(@isempty, C);
% update
if any(I),  C(I) = [];  iscol(I) = [];  ind(I) = [];  nC = numel(C);  end

%%
% check
if nC==1  % nothing more to do
    if any(iscol),  C(iscol) = cellfun(@transpose, C(iscol), 'Uniform', false);  end
    if isarray,  C = char(C);  end
    return;
end

%%
% sort the strings according to their (increasing) length
[~,I] = sort(cellfun(@length,C));  
C = C(I);  ind = ind(I);  iscol = iscol(I);

%%
% build the cell of inverted strings
if traversal<=0
    flipC = cellfun(@(c) fliplr(c), C, 'Uniform', false);
end

if traversal==1
    I = cellfun(@(c) strcompare(c, C, length(c)), C, 'Uniform', false);
elseif traversal==-1
    I = cellfun(@(flipc) strcompare(flipc, flipC, length(flipc)), ...
        flipC, 'Uniform', false);
else
    I = cellfun(@(c,flipc) ...
        strcompare(c, C, length(c)) | strcompare(flipc, C, length(c)) | ...
        strcompare(c, flipC, length(c)) | strcompare(flipc, flipC, length(c)), ...
        C, flipC, 'Uniform', false);
    % note that, in its first occurence, C is used as a constant (and,
    % similarly, flipC), so that every single string of C (second occurrence)
    % is compared with every other strings
end
I = reshape(transpose(cell2mat(I)),[nC nC]); % square matrix 

%%
% retrieve the lower triangular part of the diagonal; we do not want the
% diagonal as it is always true (the substring is compared with itself) and
% we do not want the upper triangular part as it means that two substrings
% are equal (considering that the strings were initially ordered according
% to their length) and one of them only (but not both) should be discarded
I = tril(I,-1) ;

%%
% find the strings that are substrings of longer strings already in C: they
% are repeated entries given by the indices of the non null columns of I
[~,I] = find(I); I = unique(I);

%%
% clean up: get rid of those substrings and reconvert back to the original
% format
if any(I),  C(I) = []; iscol(I) = [];  ind(I) = [];  end
if any(iscol),  C(iscol) = cellfun(@transpose, C(iscol), 'Uniform', false);  end
I = ind;

%%
% possibly reconvert back to a matrix
if isarray,  C = char(C);  end

end % end of cellstrsubtrim
