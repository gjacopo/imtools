%% LOCALGLOV2D_BASE - Base function for LOCALGLOV2D.
%
%% Syntax
%     O = LOCALGLOV2D_BASE(I, feat, res, win, wei, n, sig, mask );
%
%% See also
% Related:
% <LOCALGLCM2D_BASE.html |LOCALGLCM2D_BASE|>,
% <LOCALGLSDV2D_BASE.html |LOCALGLSDV2D_BASE|>.
% Called:
% <../../kernel/html/NEIPOSKERNEL.html |NEIPOSKERNEL|>,
% <../../kernel/html/EUCLIDKERNEL.html |EUCLIDKERNEL|>,
% <matlab:webpub(whichpath('ACCUMARRAY')) |ACCUMARRAY|>.

%% Function implementation
function O = localglov2d_base(I, feat, res, win, wei, n, sig, mask )

warning('localglov2d_base:warning', ...
    'currently computing the entropy feature only');

%% 
% parsing parameters

if ~isnumeric(I)
    disp('a matrix is required in input'); return;
end
[M,N,C] = size(I);

%% 
% checking: ensure that the window size is odd
if mod(win,2)==0,  win = win + 1;  end
ws = floor(win / 2); % half size window

%%
% dealing with multichannel images

O = zeros(M,N);
if C>1
    for c=1:C
        O(:,:,c) =  localglov2d_base(I(:,:,c), feat, res, win, wei, n, sig, mask );
    end
end

%% 
% initializing

if strcmp(n,'global');
    Imin = min(I(:));
    Idelt = max(I(:)) - Imin;
end

%%
% pad the input image
A = padarray(I, [ws ws],'replicate','both'); 
[X,Y,~] = size(A);                                                     

%%
% create the (temp) output image
O = zeros(X,Y);

%%
% indexes
pixindex = reshape(1:X*Y,X,Y);
pixin = reshape(pixindex(ws+1:ws+M,ws+1:ws+N),1,M*N);

%%
% index of neighbour window of analysis 
indI = neiposkernel(ws, X);
indI = indI(:);

%%
% initialize the weights
switch wei
    case 'ave'
        W = fspecial('average',win);
    case 'gaus'
        W = fspecial('gaussian',win,sig);
    case 'inv'
        W = euclidkernel([win win], 1, false, true);
        W(ws+1,ws+1) = 1;
        if sig~=1, W = W.^sig; end
end
W = W(:);
W = W / sum(W);

%%
% main computation
W = [0; W; 0];

for in=pixin

    % extract current values in window
    inI = A(indI + in);    
    if strcmp(n,'local');
        Imin = min(inI(:));
        Idelt = max(inI(:)) - Imin;
    end
    % quantize
    inI = ceil((res-1) * (inI - Imin) / Idelt + 0.5); % in range [1,nbin]
    
    %%
    % fast method
    % ensure that the output vector is of size nbin, and not shorter, by
    % padding the greylevels 1 and res with null probability (see W
    % padding above), so that accumarray take them into account in the 1D
    % vector construction
    inI = [1; inI; res];                                               %#ok
    % accumulate using the weights in W to create the 1D vector of size
    % nbin
    pI = accumarray(inI,W);
    % compute the feature
    O(in) =  entropy2(pI(:));    
    
end

O = reshape(O(pixin),M,N);

end % end of localglov2d_base


%% Subfunctions

%--------------------------------------------------------------------------
function E = entropy2(pI)                                              
pI(pI==0) = 1;
E = - sum(pI .* log2(pI));
end % end of entropy2


%--------------------------------------------------------------------------
function H = homogeneity(pI)                                           %#ok
H = sum(pI ./ I); 
end % end of homogeneity
