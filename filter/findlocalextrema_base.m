%% FINDLOCALEXTREMA_BASE - Base function for FINDLOCALEXTREMA.
%
%% Syntax
%     M = FINDLOCALEXTREMA_BASE(I, ord, win, graph, corr);
%
%% See also
% Related:
% <FINDLOCALEXTREMA.html |FINDLOCALEXTREMA|>,
% <FINDLOCALMAX_BASE.html |FINDLOCALMAX_BASE|>.
% Called: 
% <matlab:webpub(whichpath('ORDFILT2')) |ORDFILT2|>,
% <matlab:webpub(whichpath('STREL/GETNHOOD')) |STREL/GETNHOOD|>,
% <matlab:webpub(whichpath('MAX')) |MAX|>,
% <matlab:webpub(whichpath('SUM')) |SUM|>.

%% Function implementation
function M = findlocalextrema_base(I, ord, win, graph, corr)

% we allow variable number of inputs
if nargin<4,    graph = 8;
    if nargin<5,   corr = 'inone'; 
    end
end

%% 
% setting internal variables
[X Y C] = size(I);

% should be odd
win = round( (win-1)/2 )*2 + 1;
padnum = (win-1)/2;
A = padarray(I,[padnum padnum],'replicate','both');

%if strcmp(method,'ordfilt2')
% B=ORDFILT2(A,ORDER,DOMAIN)
if graph ==4
    v2 = getnhood(strel('diamond',win));
elseif graph==8
    v2 = true(win);
end
%vs2 = getnhood(strel('diamond',p.win));
vs2 = v2;
vs2(padnum+1,padnum+1) = 0;

switch ord
    case 'max'
        v1 = length(find(v2));
        vs1 = length(find(vs2));
    case 'min'
        v1 = 1;
        vs1 = 2;
end

% elseif strcmp(method,'nlfilter') % !!!slow method!!!
%     % B = NLFILTER(A,[M N],FUN)
%     v1 = [p.win p.win];
%     switch p.ord
%         case 'max'
%             v2 = @(x) max(x(:));
%         case 'min'
%             v2 = @(x) min(x(:));
%     end
% 
% end

%% 
% processing channel by channel

W = zeros(size(I));
for c = 1:C    
    % ie, a pixel of a multichannel image is an extrema pixel if and only
    % if it is a extrema in each individual channel.
    M = ordfilt2(A(:,:,c), v1, v2);
    Ms = ordfilt2(A(:,:,c), vs1, vs2);
    M = A(:,:,c)==M & A(:,:,c)~=Ms;
    W(:,:,c) = M(1+padnum:X+padnum,1+padnum:Y+padnum);
end

%%
% combine

if C>1 
    switch corr
        case 'inall'
            M = sum(W,3) == C;
            % or 
            % M = max(W,[],3); M = sum(M,3)>=1;
        case 'inone'
            M = max(W,[],3);
        case 'forall'
            M = W;
    end
else
    M = W;
end

end % end of findlocalextrema_base
