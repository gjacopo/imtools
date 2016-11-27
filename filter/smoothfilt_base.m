%% SMOOTHFILT_BASE - Base function for SMOOTHFILT.
%
%% Syntax
%     S = smoothfilt_base(I, rho, sm, hsize, samp, thez, sigt, theta);
%
%% See also 
% Related:
% <SMOOTHFILT.html |SMOOTHFILT|>,
% Called:
% <../../filter/html/CONVOLUTION_BASE.html |CONVOLUTION_BASE|>,
% <../../kernel/html/HOURGLASSKERNEL.html |HOURGLASSKERNEL|>,
% <../../kernel/html/GAUSSKERNEL.html |GAUSSKERNEL|>,
% <matlab:webpub(whichpath('IMGAUSSIAN')) |IMGAUSSIAN|>,
% <matlab:webpub(whichpath('FSPECIAL')) |FSPECIAL|>,
% <matlab:webpub(whichpath('IMFILTER')) |IMFILTER|>,
% <matlab:webpub(whichpath('CONV2')) |CONV2|>.

%% Function implementation
function S = ...
    smoothfilt_base(I, rho, sm, hsize, samp, thez, sigt, theta)

%% 
% checking parameters and setting variable
[X Y C] = size(I);                                                 

% original unsampled sizes
sX = X / samp; sY = Y / samp;

if isempty(hsize)
    
    hsize = ceil(6 * samp * rho);
    hsize = hsize + ~mod(hsize,2); % make it odd
    % switch sm
    %     case {'conv','convolution'}
    %         hsize = 1 + 2*round( samp * rho * 1.2 );
    %         hsize = min(min(hsize, 1+2*round(X/2)), 1+2*round(Y/2));
    %     case {'ani','hourglass'}
    %         hsize = max(1 + 2*round(rho),7);
    %     otherwise
    %         hsize = ceil(6 * samp * rho);
    % end
end

%%
% dealing with multispectral images
if C>1
    S = zeros(X,Y,C);
    for c=1:C
        S(:,:,c) = ...
            smoothfilt_base(I(:,:,c), rho, sm, hsize, samp, thez, sigt, theta);
     end
    return;
end

%% 
% linear or nonlinear spatial smoothing
    
% 4 different techniques for smoothing
switch sm
    
    % perform classical isotropic smoothing using a Gaussian filter
    case {'fast','imgaussian'}
        S = imgaussian( I, samp * rho, hsize );
        
    case {'conv','convolution'}
        h = gausskernel([hsize hsize], samp*rho/sqrt(sX*sY), [sX sY]);
        % h = gausskernel([hsize hsize], samp*rho, 1);
        S = convolution_base(I, h, 'sym');
        
    case {'matlab','conv2'}
        % S = imfilter(I, gaussian);
        % use separability of Gaussian filters for faster implementation
        gaussian = fspecial('gaussian',[1,hsize],rho); % 1D filter
        % S = imfilter(imfilter(I,gaussian','same','replicate'),gaussian,'same','replicate')
        S = conv2(gaussian,gaussian,I,'same');
        
        % perform anisotropic non-linear hour-glass smoothing
    case {'ani','glass'} 
        % preallocate
        S = zeros(X,Y,thez);
        % smoothing along the long edge direction
        vtheta = linspace(0,pi,thez+1); vtheta(end) = [];
        for i=1:thez
            h = hourglasskernel( samp*[hsize hsize], ...
                samp*rho/sqrt(sX*sY), ...
                sigt, [sX sY], vtheta(i) );
            % figure, imagesc(h), colormap gray, title(num2str(vtheta(i)))
            S(:,:,i) = convolution_base(I, h, 'sym');
        end
        % edge orientation estimation in the upper quadrant (the
        % hourglass filters are symmetric)
        theta = repmat( theta, [1 1 thez] );
        vtheta = repmat( reshape(vtheta(:),[1 1 thez]), [X Y 1] );
        % compute the min among the theta difference
        [~,M] = min( abs(theta-vtheta),[],3 );
        % select correct location
        M = reshape( (1:(X*Y))' + (M(:)-1)*(X*Y), X, Y );
        %M=reshape( (1:(X*Y))', X, Y );
        S = S(M);
        
end

end % end of smoothfilt_base
