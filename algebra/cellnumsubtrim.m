%% CELLNUMSUBTRIM - Trim a set of numerical sequences from its subsequences. 
%
%% Description
% Discard all the elements of a set (cell) of numerical sequences that match
% subsequences of longer elements in the same set.
%
%% Syntax
%       [C, I] = CELLNUMSUBTRIM(S);
%       [C, I] = CELLNUMSUBTRIM(S, traversal);
%
%% Inputs
% *|S|* : cell of numerical sequences from which redundant elements matching
%     subsequences of longer sequences are removed.
%
% *|traversal|* : (optional) flag setting the rule for defining subsequences;
%     a numerical sequence in |S| is regarded as a subsequence if it matches
%     a subsequence of a longer numerical sequence of |S|:          
%  
% * when reading both of them *from left to right* in the case |traversal=1|,          
% * when reading both of them *from right to left* in the case |traversal=-1|,  
% * *indepently of the order these strings are read* in the case |traversal=0|;     
%          
% see also |CELLSTRSUBTRIM| for string comparison; default: |traversal=0|.    
%
%% Outputs
% *|C|* : output cell of numerical sequences where the subsequences of |S|
%     have been removed following the rule defined by |traversal|; note
%     that |C| is ordered according to the length of the sequences.
%     
% *|I|* : indices (referring to the input |S|) of the numerical sequences 
%     present in the output.
%     
%% Example
%
%% See also
% Related:
% <../../misc/html/CELLSTRSUBTRIM.html |CELLSTRSUBTRIM|>.
% Called:
% <../../algebra/html/MAT2RC.html |MAT2RC|>.
% <matlab:webpub(whichpath('TRIL')) |TRIL|>,
% <matlab:webpub(whichpath('STRNCMP')) |STRNCMP|>,
% <matlab:webpub(whichpath('ISSPACE')) |ISSPACE|>,
% <matlab:webpub(whichpath('CELLSTR')) |CELLSTR|>,
% <matlab:webpub(whichpath('CELL2MAT')) |CELL2MAT|>,
% <matlab:webpub(whichpath('CELLFUN')) |CELLFUN|>,
% <matlab:webpub(whichpath('RESHAPE')) |RESHAPE|>.

%% Function implementation
%--------------------------------------------------------------------------
function [C, I] = cellnumsubtrim(C, traversal)

%%
% check/set variables

if nargin<2,  traversal = 0;  end

if ~iscell(C)
    error('cellnumsubtrim:errorinput', ...
        'input variable ''C'' must be a cell of vectors');
    
elseif ~isscalar(traversal) || ~ismember(traversal,[0 1 -1])
    error('cellnumsubtrim:errorinput', ...
        'input variable ''traversal'' takes values in {-1,0,1}');
end

nC = numel(C); 
ind = 1:nC;

%%
% in general, get rid of empty sequences
I = cellfun(@isempty, C);
% update
if any(I),  C(I) = [];  ind(I) = [];  nC = numel(C);  end

%%
% check if one sequence only: then nothing to do
if nC==1,  return;  end

%%
% sort the sequences according to their (increasing) length
[~,I] = sort(cellfun(@length,C));  
C = C(I);  ind = ind(I);

%%
% now we use a trick to get rid of redundant sequences: we transform the
% numerical sequences into strings and we compare the strings using the
% standard |STRNCMP| function, likewise we do in |CELLSTRSUBTRIM|
A = cellfun(@(c) num2str(mat2rc(c,'r')), C, 'Uniform', false);
% note that we also ensure that the (string) sequences are row vectors...

%%
% build the cell of inverted strings
if traversal<=0
    flipA = cellfun(@(c) num2str(fliplr(mat2rc(c,'r'))), C, 'Uniform', false);
end

%%
% clean up a bit: first deblank (remove leading and trailing whitespace
% characters from string, not optional), then possibly'despace' (trim all
% internal whitespaces)
A = cellfun(@strtrim, A, 'Uniform', false);
if traversal<=0,  flipA = cellfun(@strtrim, flipA, 'Uniform', false);  end

%%
% check the presence of substrings
if traversal==1
    I = cellfun(@(a) strncmp(a, A, length(a)), A, 'Uniform', false);
elseif traversal==-1
    I = cellfun(@(flipa) strncmp(flipa, flipA, length(flipa)), ...
        flipA, 'Uniform', false);
else
    I = cellfun(@(a,flipa) ...
        strncmp(a, A, length(a)) | strncmp(flipa, A, length(a)) | ...
        strncmp(a, flipA, length(a)) | strncmp(flipa, flipA, length(a)), ...
        A, flipA, 'Uniform', false);
    % note that, in its first occurence, A is used as a constant (and,
    % similarly, flipA), so that every single string of A (second occurrence)
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
[~,I] = find(I);  I = unique(I);

%%
% clean up: get rid of those substrings and reconvert back to the original
% format
if any(I),  C(I) = [];  ind(I) = [];  end
I = ind;

end % end of cellnumsubtrim
