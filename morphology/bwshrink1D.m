%% BWSHRINK1D - Shrinking of 1D binary signals.
%
%% Description
% Shrink a 1D binary signal by reducing sequence of points with true value 
% to isolated true points. Similarly to |BWMORPH| (which it uses when the Image
% toolbox is available), it removes points so that a true sequence shrinks  
% to a point halfway of it.
%
%% Syntax
%       S = BWSHRINK1D(X);
%
%% Input
% *|X|* : 1D binary signal; if not binary, it is converted using |LOGICAL|.
%
%% Output
% *|S|* : output shrinked signal with same size as the input |X|.
%
%% See also
% Called:
% <matlab:web(whichpath('bwmorph')) |BWMORPH|>,
% <matlab:web(whichpath('diff')) |DIFF|>,
% <matlab:web(whichpath('find')) |FIND|>.

%% Function implementation
function S = bwshrink1D(X)

%% 
% check the inputs
if ~islogical(X),  X = logical(X);  end

if nb_dims(X)~=1, 
    error('bwshrink1D:errorinput',...
        '1D signals required in input - use BWMORPH when available');
end

%%
% in the case the Image toolbox is available, simply call the morphological
% function BWMORPH
if ~isempty(ver('images'))
    S = bwmorph(X, 'shrink', Inf);
    return
end 
% else...

%%
% initialize the output 
S = false(size(X));

%%
% characterize the logical changes in sequences
x1 = [X(1); diff(X(:))];
x2 = [-diff(fliplr(X(:))); X(end)];

%%
% extract those isolated 'true' points that are shrinked on themselves
isolated = find(x1==1 & x2==1);
S(isolated) = true;

%%
% extract the starting and ending points of a 'true' sequence
x1(isolated) = 0;  x1 = find(x1==1);
x2(isolated) = 0;  x2 = find(x2==1);

%%
% shrink the sequence to its central point
S(ceil((x1+x2)/2)) = true; % with CEIL, returns the same output as BWMORPH
% S(floor((x1+x2)/2)) = true;

%%
% final output
S = reshape(S,size(X));

end % end of bwshrink1d

