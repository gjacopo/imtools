%% WORDFILT2_BASE - Base function for WORDFILT2.
% 
%% Syntax
%     WOI = wordfilt2(I, wk);
%     WOI = wordfilt2(I, wk, order)
%
%% See also
% Related:
% <WORDFILT2.html |WORDFILT2|>,
% <WMEDFILT2_BASE.html |WMEDFILT2_BASE|>.
% Called:
% <matlab:webpub(whichpath('SORT')) |SORT|>,
% <matlab:webpub(whichpath('MEDIAN')) |MEDIAN|>,
% <matlab:webpub(whichpath('MIN')) |MIN|>,
% <matlab:webpub(whichpath('MAX')) |MAX|>,
% <matlab:webpub(whichpath('MODE')) |MODE|>.

%% Function implementation
function WOI = wordfilt2_base(I,wk,order)

[X Y C] = size(I);
[x y c] = size(wk);

cpos = find(wk<0,1,'first');
if isempty(cpos),  cpos = (x+1+x*(y-1))/2; end;

[i,j] = ind2sub([x y],cpos); % position of the center

%%
% dealing with multispectral matrices

if C>1
    WOI = zeros(size(I));
    for ic=1:C
        if c~=1; icc=ic; else icc=1; end
        WOI(:,:,ic) = wordfilt2_base(I(:,:,ic), wk(:,:,icc), order);
    end
    return
end

%% 
% main processing

% pad the image
A = padarray(I,[i-1 j-1],'symmetric','pre');
A = padarray(A,[x-i y-j],'symmetric','post');
Xpad = X+x-1; % size(A,1)
Ypad = Y+y-1; % size(A,2)

% index of the pixels inside the padded image
apix = reshape(1:Xpad*Ypad, [Xpad Ypad]);
ipix = apix(i:X+i-1,j:j+Y-1);

% index of the kernel
ik = repmat((1:y)-j,[x 1])*Xpad + repmat(((1:x)-i)',[1 y]);

%%
% create the 3D matrice with weightened values: the k-th pixel value in the 
% neighbourhood specified by wk is replicated along the 3rd dimension wk(k)
% times
WI = [];
for k=1:x*y
    if wk(k)
        % replicate wk(k) times along the 3rd dimension
        WI = cat(3,WI,repmat(A(ipix+ik(k)),[1 1 abs(wk(k))]));
    end
end

%%
% compute the order statistics
if isscalar(order)
    WI = sort(WI,3);
    WOI = WI(:,:,order);
    
else % ischar(order)
    switch order
        case 'med'
            WOI = median(WI,3);
                    
        case 'mod'
            WOI = mode(WI,3);
    end
end

end % end of wordfilt2_base
