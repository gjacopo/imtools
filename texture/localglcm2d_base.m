%% LOCALGLCM2D_BASE - Base function for LOCALGLCM2D.
%
%% Syntax
%     O = LOCALGLCM2D_BASE(I, feat, res, dcar, win, wei, sig, mask);
%
%% See also
% Related:
% <LOCALGLOV2D_BASE.html |LOCALGLOV2D_BASE|>,
% <LOCALGLSDV2D_BASE.html |LOCALGLSDV2D_BASE|>.
% <../../kernel/html/NEIPOSKERNEL.html |NEIPOSKERNEL|>.
% Called:
% <HISTCONTRAST.html |HISTCONTRAST|>,
% <HISTVARIANCE.html |HISTVARIANCE|>,
% <HISTENERGY.html |HISTENERGY|>,
% <HISTMAXIMUM.html |HISTMAXIMUM|>,
% <HISTENTROPY2.html |HISTENTROPY2|>,
% <HISTENTROPY10.html |HISTENTROPY10|>,
% <HISTDISSIMILARITY.html |HISTDISSIMILARITY|>,
% <HISTIDIFFERENCE.html |HISTIDIFFERENCE|>,
% <HISTCORRELATION.html |HISTCORRELATION|>,
% <HISTHOMOGENEITY.html |HISTHOMOGENEITY|>,
% <matlab:webpub(whichpath('ACCUMARRAY')) |ACCUMARRAY|>.
 
%% Function implementation
function O = localglcm2d_base(I, feat, res, dcar, win, wei, sig, mask)

%%
% checking parameters and setting internal variables

[M,N,C] = size(I);

%%
% ensure that the window size is odd
if mod(win,2)==0,  win = win + 1;  end
ws = floor(win / 2); % half size window

%%
% use the Cartesian representation in the following
dnorm = sqrt(dcar(:,1).*dcar(:,1) + dcar(:,2).*dcar(:,2));         
nd = length(dnorm); % also: size(dcar,1);

if strcmp(feat,'all'),    nfeat = 8;
else                      nfeat = 1;
end

%%
% dealing with multichannel images

if C>1
    O = zeros(M,N,C,nfeat,nd);
    for c=1:C
        O(:,:,c,:,:) = localglcm2d_base(I(:,:,c), feat, res, dcar, win, wei, sig, mask );
    end
    return;
end

%%
% computation

% pad the input image
pad = ws + ceil(max(dnorm));
A = padarray(I, [pad pad],'replicate','both'); 
[X,Y] = size(A);
A = A(:);

% global normalization & quantization
%if strcmp(n,'global');
Amin = min(A); Adelt = max(A) - Amin;
A = (A - Amin) / Adelt; % in range [0,1]
A = ceil((res-1) * A + 0.5);
%end

% indexes
pixindex = reshape(1:X*Y,X,Y);
pixin = reshape(pixindex(pad+1:pad+M,pad+1:pad+N),1,M*N);

% Index of the centered neighbour window of analysis 
indI = neipos(ws, X);
indI = indI(:);
% indI = -ws:ws;
% for i=1:ws
%     indI = [ (-i*X-ws):(-i*X+ws), indI, (i*X-ws):(i*X+ws) ];         %#ok
% end
% indI = indI';
% indI = indI(:);

% Index of the displaced window(s)
indJ = zeros(length(indI),nd);
for d=1:nd
    indt = dcar(d,2)*X + dcar(d,1);
    indJ(:,d) = indI + indt;
end

%%
% initialize the set of estimated features and their number

switch feat
    case 'all'
        histfeat = @histfeatures;
    case 'var',           histfeat = @histvariance;
    case 'con',           histfeat = @histcontrast;
    case 'ene',           histfeat = @histenergy;
    case 'max',           histfeat = @histmaximum;
    case 'dis',           histfeat = @histdissimilarity;
    case 'inv',           histfeat = @histidifference;
    case 'mean',          histfeat = @histmean;
    case {'ent2','ent'},  histfeat = @histentropy2;
    case 'ent10',         histfeat = @histentropy10;
    case 'cor',           histfeat = @histcorrelation;
    case 'hom',           histfeat = @histhomogeneity;
end

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
[I,J] = meshgrid(1:res);
% I: index of column [1:res] repeated on res rows
%     I = [ 1 2 3...res; 1 2 3...res; ... ; 1 2 3...res ]
% J: index of row [1:res] repeated on res columns
%     J = [ 1 1 1...1; 2 2 2...2; ... ; res res res...res ]
I = I(:); J = J(:);
W = [0; W; 0];

% create the (temp) output image
O =  zeros(X*Y,nfeat,nd);
  
%% 
% loop over the image
for in=pixin

    %%
    % extract current values in windows
    inI = A(indI + in);
    % normalize locally if required
    % if strcmp(n,'local');
    %    mmin = min(inI(:));
    %    inI = (inI - mmin) / (max(inI) - mmin); % in range [0,1]
    % end
    % inI = ceil((res-1) * inI + 0.5); % in range [1,res]
    inI = [1; inI; res];                                               %#ok
    
    %%
    % loop over the displacement vectors
    for d=1:nd
        inJ = A(indJ(:,d) + in);
        %  inJ(:,d) = ceil((res-1) * inJ(:,d) + 0.5);
        % if strcmp(n,'local');
        %    mmin = min(inJ(:,d));
        %    inJ(:,d) = (inJ(:,d) - mmin) / (max(inJ(:,d)) - mmin);
        % end
        
        %%
        % fast method
        % ensure that the output matrix is of size (res x res) by padding the
        % greylevels pairs (1,1) and (res,res) with null probability (see W
        % padding above), so that accumarray take them into account in the 2D
        % array construction
        inJ = [1; inJ; res];                                           %#ok 
        % accumulate using the weights in W to create the 2D matrix of size
        % (res x res)
        pIJ = accumarray([inI inJ],W);
        %      pIJ0 = pIJ(:)==0;
        %      pIJ(pIJ0) = [];
        %      II = I; II(pIJ0) = [];
        %      JJ = J; JJ(pIJ0) = [];
        % compute the feature
        O(in,:,d) = histfeat(pIJ(:),I,J);
         
        % elegant method...but slower
        %      pIJ = accumarray([inI inJ],W);
        %      IJ = unique([inI inJ],'rows');
        %      pIJ(pIJ==0) = [];
        %      O(in) =  histvariance(pIJ,IJ(:,1),IJ(:,2));
    end

end

O = reshape(O(pixin,:,:),M,N,nfeat,nd);

end % end of localglcm2d_base
