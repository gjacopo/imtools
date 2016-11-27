%% POSINSERT - Insert numeric values in a vector or matrix at given positions.
% 
%% Syntax
%       V = posinsert(v, pos);
%       V = posinsert(v, pos, n, flag);
%
%% Inputs
% *|v|* : input vector (or matrix) to pad by inserting values at given  
%     positions.
%
% *|pos|* : index positions (preferably ordered) where to insert the values
%     of |n| in |v|; if empty, the input variable |v| is returned in |V|;
%     see also |flag| below.
%     
% *|n|* : (optional) value(s) to insert in |v|; it must be of size |numel(pos)|
%     or smaller; in this latter case, it is padded so to fit |pos|; in
%     particular, if a scalar is passed, it will be inserted everywhere as
%     specified by |pos|; default: |n=0|.
% 
% *|flag|* : (optional) string specifying if the insertion position is 
%     relative to the input variable or the output variable:
%
% * in the case |'i'|, the |i|-th value of |n| will appear in the output vector
%          between the |i-1|-th and the |i|-th values of the input vector; it is
%          in fact inserted at the position |pos(i)+i-1|;
% * in the case |'o'|, the |i|-th value of |n| appears (as expected...) at the
%          position |pos(i)| in the output variable;
%
% default: |flag='i'|.     
%     
%% Output
% *|V|* : padded variable where values of |n| have been inserted in |v|
%     according to the positions specified in |pos|.
%     
%% Example
%  a=1:10;  pos=[4 7 8];
%  posinsert(a,pos,0,'i') 
%  posinsert(a,pos,[20 30 40],'o') 
%  
%
%% Remark
% Calling |posinsert(v, pos, n, 'i')| is in fact strictly equivalent to calling
% |posinsert(v, pos+(0:length(pos)-1), n, 'o')|.  

%% Function implementation
%--------------------------------------------------------------------------
function V = posinsert(v, pos, n, flag)

%%
% check variables

% dummy case
if isempty(pos), 
    V = v; 
    return;
end

if nargin<4,  flag = 'i'; 
    if nargin<3 || isempty(n),  n = 0;  end
end

if ~any(strcmpi(flag,{'i','o'}))
    error('posinsert:errorinput', 'unknown flag - see help')
end

if numel(n)>numel(pos) || ...
        (strcmp(flag,'i') && max(pos)>numel(v)) || ...
        (strcmp(flag,'o') && max(pos)>numel(v)+numel(pos))
    error('posinsert:errorinput', ...
        'check the size compatibility of pos with the input vector')
end

%%
% set default values

if numel(n)<numel(pos)
    n = padarray(n, numel(pos)-numel(n), n(end), 'post'); 
    % in, particular if one value is passed only, this value is inserted
    % everywhere at pos positions
end

%%
% positions of appearance should be ordered in increasing order! so reorder,
% just in case...
[pos, I] = sort(pos(:));
n = n(I);

%%
% deal with the case 'position relative to the input'
if strcmp(flag,'i'),   pos = pos(:) + (0:length(pos)-1)';  end

%%
% main processing: insert
tf = false(1,numel(v)+numel(n));

V = zeros(size(v));
tf(pos) = true;

V(tf) = n;
V(~tf) = v;
end % end of posinsert
