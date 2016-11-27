%% WMEDFILT2_BASE - Base function for WMEDFILT2.
% 
%% Syntax
%     WMI = WMEDFILT2_BASE(I, wk);
%
%% See also
% Related:
% <WMEDFILT2.html |WMEDFILT2|>,
% <WORDFILT2_BASE.html |WORDFILT2_BASE|>.
% Called:
% <matlab:webpub(whichpath('MEDIAN')) |MEDIAN|>.

%% Function implementation
function WMI = wmedfilt2_base(I,wk)

[X Y C] = size(I);
[x y c] = size(wk);

cpos = find(wk<0,1,'first');
if isempty(cpos),  cpos = (x+1+x*(y-1))/2; end;

[i,j] = ind2sub([x y],cpos);

%% 
% dealing with multispectral matrices
WMI = zeros(size(I));

% loop over C components
if C>1
    for ic=1:C
        if c~=1; icc=ic; else icc=1; end
        WMI(:,:,ic) = wmedfilt2_base(I(:,:,ic), wk(:,:,icc));
    end
    return
end

%%
% main processing

% pad the image
A = padarray(I,[i-1 j-1],'symmetric','pre');
A = padarray(A,[x-i y-j],'symmetric','post');
Xpad = X+i; Ypad = Y+j;

% index of the pixels inside the padded image
apix = reshape(1:Xpad*Ypad,Xpad,Ypad);
ipix = apix(i:X+i-1,j:j+Y-1);

% index of the kernel
ik = repmat((1:y)-j,3,1)*Ypad +repmat(((1:x)-i)',1,3);

WMI = [];
for k=1:x*y
    if wk(k)
        % replicate along the third dimension
        WMI = cat(3,WMI,repmat(A(ipix+ik(k)),[1 1 abs(wk(k))]));
    end
end

WMI = median(WMI,3);

end % end of wmedfilt2_base
