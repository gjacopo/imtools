%% BLURMAP_BASE - Base function for BLURMAP.
%
%% Syntax
%     [Iblur, qmap] = BLURMAP_BASE(I, blur, qstep);
% 
%% See also
% Related:
% <BLURMAP.html |BLURMAP|>.
% Called: 
% <matlab:webpub(whichpath('FSPECIAL')) |FSPECIAL|>,
% <matlab:webpub(whichpath('IMFILTER')) |IMFILTER|>.

%% Function implementation
function [Iblur,varargout] = blurmap_base(I, blur, qstep)

[X,Y,C] = size(I);

%%
% dealing with multispectral images

if C>1
    Iblur = zeros(X, Y, C);
    if nargout==2,  varargout{1} = zeros(X, Y, C);    end
    for c=1:C
        [Iblur(:,:,c) tmp] = blurmap_base(I(:,:,c), blur, qstep);
    end
    if nargout==2,  varargout{1}(:,:,c) = tmp;    end
    return;
end

%%
% quantize the bluring map by adjusting the range of blur parameters
blur(blur<qstep) = qstep;
quantblur = qstep * round(blur/qstep);

%%
% get the list of possible blur parameters
lblur = unique(quantblur(:));
nblur = length(lblur); % no of bluring parameters

%%
% this function corrects for boundary conditions by  reflecting the image
% at the boundaries before performing the isotropic gaussian filter.
% This eliminates the "image darkening" at the edges.
Iblur = zeros(size(I(:,:,1)));

for i=1:nblur
    sigma = lblur(i);
    G = fspecial('gaussian',fix(6*sigma),sigma);
    Ifilt = imfilter(I,G,'symmetric','same');
    imap = quantblur == sigma;
    Iblur(imap) = Ifilt(imap);
end

if nargout==2
    varargout{1} = quantblur;
end

end % end of blurmap_base
