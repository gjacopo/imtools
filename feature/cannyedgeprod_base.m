%% CANNYEDGEPROD_BASE - Base function for CANNYEDGEPROD.
%
%% Syntax
%   [edgemap, mag, or] = CANNYEDGEPROD_BASE(I, sigma2, sigma1, der, c, reduce );
%
%% See also
% Related:
% <CANNYEDGEPROD.html |CANNYEDGEPROD|>,
% <CANNYEDGE_BASE.html |CANNYEDGE_BASE|>.
% Called:
% <../../derive/html/GRDSMOOTH_BASE.html |GRDSMOOTH_BASE|>,
% <CANNYEDGEMAP_BASE.html |CANNYEDGEMAP_BASE|>. 

%% Function implementation
function [edgemap, mag, or] = cannyedgeprod_base(I, sig2, sig1, der, c, reduce)

[X,Y,C] = size(I);

%%
% ensure to combine two distinct scales
if sig1==sig2,  sig2 = sig1 + 0.05;  end;
if sig1>sig2,  tmp = sig1;  sig1 = sig2;  sig2 = tmp;  end

%%
% deal with the special case of multispectral images
if C>1
    if islogical(reduce) && ~reduce % process channel by channel
        edgemap = false(size(I));
        mag = zeros(size(I));        or = zeros(size(I));
        for c=1:C
            [edgemap(:,:,c), mag(:,:,c), or(:,:,c)] = ...
                cannyedgeprod_base(I(:,:,c), sig2, sig1, der, c, reduce);
        end
        return;
        
    elseif islogical(reduce) || ...
            (ischar(reduce) && any(strcmpi(reduce,{'isum','igray','imax'})))
        % at that point, if reduce is logical, it is necessarly true
        
        if C==3 && strcmpi(reduce,'igray')
            [I,d] = rescale(I,0,1);
            I = rgb2gray(I);
            
        elseif islogical(reduce) || strcmpi(reduce,'isum')
            if C==3
                I = 0.29*I(:,:,1) + 0.59*I(:,:,2) + 0.11*I(:,:,3);
                [I,d] = rescale(I,0,1);
            else
                [I,d] = rescale(I,0,1);
                I = sum(I, 3);
            end
            
        elseif strcmpi(reduce,'imax')
            [I,d] = rescale(I,0,1);
            I = max(I, [], 3);
        end
        
        %%
        % also rescale the adjustment parameter accordingly
        c = c / d;
        
    end
end

%%
% compute the 2D filter responses at finest scale: |sig1|
[gx1,gy1] = grdsmooth_base(I, sig1, der, [], 'xy');

%%
% compute responses at coarsest scale: |sig2|
[gx2,gy2] = grdsmooth_base(I, sig2, der, [], 'xy');

if C>1 && ischar(reduce) && strcmpi(reduce,'gmax');
    gx2 = max(gx2, [], 3); gy2 = max(gy2, [], 3);
    gx1 = max(gx1, [], 3); gy1 = max(gy1, [], 3);
end

%%
% compute the correlation coefficients between the directional filter
% responses according to Sec.3.5
rhox = corr(gx1(:),gx2(:));  
rhoy = corr(gy1(:),gy2(:));

rho = sqrt(2^3 * sig1^3 *sig2^3 / (sig1^2+sig2^2)^3);              %#ok

%%
% derive the threshold according to Eq.(3.12) and Sec.3.6
kappax = sqrt(1+2*rhox*rhox) * sig1 * sig2; 
kappay = sqrt(1+2*rhoy*rhoy) * sig1 * sig2;
thresp = sqrt(c * (kappax + kappay)); 

%%
% define the scale product responses
gx2 = gx1 .* gx2;
gx2(gx2<0) = 0;
gy2 = gy1 .* gy2;
gy2(gy2 < 0) = 0;

%%
% define the combined gradient magnitude according to Eq.(3.13)
mag = sqrt(gx2 + gy2);

%%
% threshold to remove noise
magmax = max(mag(:));
if magmax>0
    mag = mag / magmax;   % normalize
    thresp = thresp / magmax;
end
mag(mag<=thresp) = 0;

%%
% define the combined gradient orientation according to Eq.(3.13); the
% orientation information is recovered from the scale product responses and
% the finest scale |sig1| (through the signs of |gy1| and |gx1|)
or = atan2(sign(gy1) .* sqrt(gy2), sign(gx1) .* sqrt(gx2)); 
%or = atan2(-sqrt(gx1) .* sign(gx1),sqrt(gy2) .* sign(gy1)); 

%%
% finally compute the edge map
edgemap = cannyedgemap_base(gx2, gy2, der, mag, or, [], []);
% ratio = [1/3 0.08];
% hyst = [];

%%
% final output, special case
if C>1 && ischar(reduce) && strcmpi(reduce,'eor')
    [edgemap, J] = max(edgemap, [], 3);
    subm = reshape(1:X*Y,[X Y]) + (J-1)*X*Y;
    mag = mag(subm);
    or = or(subm);
end

end % end of cannyedgeprod_base
