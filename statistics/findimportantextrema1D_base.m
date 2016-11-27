%% FINDIMPORTANTEXTREMA1D_BASE - Base function for FINDIMPORTANTEXTREMA1D.
%
%% Syntax
%    [Extrema, Strict, Left, Right] = findimportantextrema1D_base(X, dist);
%
%% See also
% Related:
% <FINDZEROEXTREMA1D_BASE.html |FINDZEROEXTREMA1D_BASE|>,
% Called:
% <FINDZEROEXTREMA1D_BASE.html |FINDZEROEXTREMA1D_BASE|>,
% <matlab:webpub(whichpath('BWCONNCOMP')) |BWCONNCOMP|>,
% <matlab:webpub(whichpath('REGIONPROPS')) |REGIONPROPS|>,
% <matlab:webpub(whichpath('LABELMATRIX')) |LABELMATRIX|>.

%% Function implementation
function [Extrema, Strict, Left, Right] = findimportantextrema1D_base(X, dist)

%%
% check/set parameters

switch dist
    % dist = @(a,b) abs(a-b) / max(abs(a),abs(b));
    case 'abs'
        dist = @(a,b) abs(a-b);
        
    case 'nabs'
        dist = @(a,b) abs(a-b) / (abs(a)+abs(b));
end

Extrema = zeros(size(X),2);

%%
% first identify the local (flat or strict) extrema

% function handle for finding (flat or strict) minima
Min = findzeroextrema1D_base(X(:), [], [], 'local3', [], [eps 0]); % [eps 0]
Extrema(Min(:,1),:) = [ones(size(Min,1),1) Min(:,2)];
% setting A = diff(CSP(1:end-1)) and B = diff(CSP(2:end)), we look for those
% edges verifying:
%    sign(A)-sign(B)<0 & min(abs([A B]),[],2)>=eps & max(abs([A B]),[],2)>0
[Strict,Left,Right] = importance(X(:), Min);

%   %----------------------------------------------------------------------
    function [strict,left,right] = importance (X, Mset)
        left = zeros(size(X));
        right = zeros(size(X));
        
        for v=unique(Mset(:,2))'
            % look at extrema with same level (ie. value)
            m = Mset(Mset(:,2)==v,1); 
            CC = bwconncomp(X>v);  L = labelmatrix(CC);
            STATS = regionprops(CC, X, 'MaxIntensity');
            maxXval = cat(1, STATS.MaxIntensity);
            
            l = L(m-1); i = l~=0; l = l(i);           
            left(m(i)) = dist(maxXval(l), v);
            left(m(~i)) = NaN;

            r = L(m+1); i = r~=0; r = r(i);
            right(m(i)) = dist(maxXval(r), v);
            right(m(~i)) = NaN;
            
        end
        strict = max([left right],[],2);
        
    end
%   %----------------------------------------------------------------------

%%
% ibid with max
Max = findzeroextrema1D_base(X(:), [], 'local3', [], [], [eps 0]); % [eps 0]
Extrema(Max(:,1),:) = [ones(size(Max,1),1) Max(:,2)];
Max(:,2) = max(X(:)) - Max(:,2);
[strict,left,right] = importance(max(X(:))-X(:), Max);
Left = Left + left;
Right = Right + right;
Strict = Strict + strict;

%%
% output
v = find(Extrema(:,1));  Extrema = [v Extrema(v,2)];
Strict = nonzeros(Strict);
Left = nonzeros(Left);
Right = nonzeros(Right);
end % end of findimportantextrema1D_base
