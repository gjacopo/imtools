%% HARRISCORNER_BASE - Base function for HARRISCORNER.
%
%% Syntax
%     [ptcorner, cornermap] = harriscorner_base(I, sig, rho, kappa, thres, radius);
%     [ptcorner, cornermap] = harriscorner_base(gx, gy, rho, kappa, thres, radius);
%     [ptcorner, cornermap] = harriscorner_base(gx2, gy2, gxy, kappa, thres, radius);
%
%% See also
% Related:
% <HARRISCORNER.html |HARRISCORNER|>,
% <SUSANCORNER_BASE.html |SUSANCORNER_BASE|>,
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>,
% <FASTCPDA_BASE.html |FASTCPDA_BASE|>,
% <FASTCORNER_BASE.html |FASTCORNER_BASE|>.
% Called: 
% <../../statistics/html/FINDLOCALMAX_BASE.html |FINDLOCALMAX_BASE|>,
% <matlab:web(whichpath('CONV2')) |CONV2|>.

%% 
% dealing with multispectral images

%% Function implementation
function [pt, varargout] = harriscorner_base(I, V1, V2, kappa, thres, radius)

[X,~,C] = size(I);                                                     

pt = cell(C,1);

if C>1
    if isscalar(V1), V1 = repmat(V1, [1 1 3]);
    elseif isempty(V1), V1 = zeros(1,1,3);   end % dummy variable
    if isscalar(V2), V2 = repmat(V2, [1 1 3]);
    elseif isempty(V2), V2 = zeros(1,1,3);   end % dummy variable
    for c=1:C
        tmp = harriscorner_base(I(:,:,c), V1(:,:,c), V2(:,:,c), kappa, radius);
        pt{c} = tmp{1};
    end
    return;
end

%%
% setting internal variables

if isscalar(V1),     sigma = V1;  
elseif isempty(V1),  sigma=1;  % also set sigma=1 in case...
end
if isscalar(V2),     rho = V2;    
elseif isempty(V2),  rho = 1;
end;

%%
% initial image filtering in the case the 2nd directional derivatives are
% not directly passed to the function

if isempty(V2) || (isscalar(V2) && V2>0)
    
    % compute the 1st-order directional derivatives
    if isscalar(V1) && V1>0
        s_D = 0.7 * sigma;
        x  = -round(3*s_D):round(3*s_D);
        dx = x .* exp(-x.*x/(2*s_D*s_D)) ./ (s_D*s_D*s_D*sqrt(2*pi));
        dy = dx';
        % compute the directional derivatives
        V1 = conv2(I, dx, 'same');
        V2 = conv2(I, dy, 'same');
    end
    
    % compute the 2nd-order directional derivatives
    % Gaussian filter
    hsize = fix(6*rho+1);
    g = fspecial('gaussian', max(1,hsize), rho);
    % compute the auto-correlation matrix from the smoothed squared image
    % derivatives: nothing else than gradient structure tensor
    I = conv2(V1.*V2, g, 'same');
    V1 = conv2(V1.^2, g, 'same');
    V2 = conv2(V2.^2, g, 'same');

end

%% 
% compute the interest point response
 
if kappa == 0
    % improved Noble/Forster measure
    cim = (V1.*V2 - I.^2)./(V1 + V2 + eps);
    
else
    % original Harris measure
    cim = (V1.*V2 - I.^2) - kappa*(V1 + V2).^2;
    % det(T) - kappa * trace(T)
end

%%
% final filtering 

% find local maxima on (radius x radius) neighbourhood
max_local = findlocalmax_base(cim, radius, 'filt');  
% set threshold 0.5% of the maximum value
thres = thres * max(max_local(:));
% find local maxima greater than threshold
 [r,c] = find(max_local>=thres);
% build interest points
pt{1} = [r,c];

if nargout == 2
    varargout{1} = false(size(I));
    varargout{1}(r+X*(c-1)) = 1;
end

end % end of harriscorner_base



