%% GRDSMOOTH_BASE - BASE function for GRDSMOOTH. 
%
%% Syntax
%      [gx, gy] = GRDSMOOTH_BASE(I, sigma, der, hsize, axis);
%      [gx, gy, mag, or] = GRDSMOOTH_BASE(I, sigma, der, hsize, axis);
%
%% Acknowledgment
% This function is a copy/paste and crop of several proposed functions 
% implemented for smoothing and differentiation following Canny's principles.
%
%% Remarks
% * See C.Luengo discussion on Gaussian filtering and Gaussian derivation:
%    http://www.cb.uu.se/~cris/blog/index.php/archives/22
%    http://www.cb.uu.se/~cris/blog/index.php/archives/150#more-150
%
% * See P.Kovesi toolbox on improvements for local gradient estimation:
%    http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/
%
%% See also  
% Related:
% <GRDSMOOTH.html |GRDSMOOTH|>.
% Called:
% <../../kernel/html/GAUSSKERNEL.html |GAUSSKERNEL|>,
% <../../filter/html/CONVOLUTION.html |CONVOLUTION|>,
% <matlab:webpub(whichpath('DERIVATIVES')) |DERIVATIVES|>,
% <matlab:webpub(whichpath('IMGAUSSIAN')) |IMGAUSSIAN|>,
% <matlab:webpub(whichpath('FSPECIAL')) |FSPECIAL|>,
% <matlab:webpub(whichpath('IMFILTER')) |IMFILTER|>,
% <matlab:webpub(whichpath('FILTER2')) |FILTER2|>,
% <matlab:webpub(whichpath('GRADIENT')) |GRADIENT|>.

%% Function implementation
function [gx,gy,varargout] = grdsmooth_base(I, sigma, der, hsize, axis)

%%
% dealing with multispectral images
C = size(I,3);

if C>1
    gx = zeros(size(I)); gy = zeros(size(I));
    for i=1:nargout-2
        varargout{i} = zeros(size(I));
    end
    for c=1:C
        [gx(:,:,c) gy(:,:,c), tmpmag, tmpor] = ...
                grdsmooth_base(I(:,:,c), sigma, der, hsize, axis);
        if nargout>=3,  varargout{1}(:,:,c) = tmpmag;  end
        if nargout==4,  varargout{2}(:,:,c) = tmpor;  end
    end
    return;
end


%% 
% handling the function

% define the string der containing the method and convert it to the 
% function handle
der = ['grdsmooth_' der];
der = str2func(der);

%%
% main computation

% estimation of the directional derivatives
[gx,gy] = der(I,sigma,hsize);                                          %#ok 

if strcmp(axis,'xy')
    % take the vector orthogonal to the output gradient to get the 'real' 
    % gradient
    tmp = gx; gx = gy; 
    % negate the second component of the vector in order to turn a left handed
    % vector (the usual result of the gradient filter, because gy runs from
    % top to bottom) into a right handed tensor
    gy = -tmp;
end

if nargout>=3
    % norm of the gradient (Combining the I and Y directional derivatives)
    mag = hypot(gx,gy);
    magmax = max(mag(:));
    if magmax>0
        mag = mag / magmax;   % normalize
    end
    varargout{1} = mag;
    
    if nargout==4
        % orientation
        if strcmp(axis,'xy'),  orien = atan2(gy, gx);           
        else                   orien = atan2(-gx, gy); 
        end  % angles -pi to + pi.
        % neg = orien<0;                   
        % varargout{2} = orien.*~neg + (orien+pi).*neg; % map angles to 0-pi.
        % orien = orien*180/pi;               % convert to degrees.
        varargout{2} = orien;
    end
end

end % end of grdsmooth_base


%% Subfunctions

%%
% |GRDSMOOTH_MATLAB| - The directional gradients' estimation using the 2D
% Gaussian filtering of the image with |IMFILTER| prior to its differentiation
% with |GRADIENT|.
%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_matlab(I, sigma, hsize)                   %#ok

if sigma>0.05
    if isempty(hsize)
        hsize = round(6*sigma)+1;   % the filter size.
    end
    
    % smooth the image out
    gaussian = fspecial('gaussian',hsize,sigma);
    Isigma = imfilter(I, gaussian);
    
else
    Isigma = I;
end

% calculate the derivatives
[gy,gx] = gradient(Isigma);

% % other approach:
% Isobel = fspecial('sobel');
% gx = -imfilter(I,Isobel);
% gy = -imfilter(I,Isobel');
end


%%
% |GRDSMOOTH_FAST| - The directional gradients' estimation using a fast
% 2D Gaussian convolution of the image with |IMGAUSSIAN| prior to its
% differentiation with |DERIVATIVES|.
%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_kroon(I, sigma, hsize)                   %#ok
[gx,gy] = grdsmooth_fast(I, sigma, hsize);
end % here for convenience (in some calls of the function GRDSMOOTH_BASE)
%--------------------------------------------------------------------------
function [ux,uy] = grdsmooth_fast(I, sigma, hsize)                   

if sigma>0.05
    if isempty(hsize)
        hsize = round(6*sigma);   % the filter size.
    end
    
    % smooth the image out
    Isigma = imgaussian(I,sigma,hsize);
else
    Isigma = I;
end

% calculate the derivatives
ux = derivatives(Isigma,'x');
uy = derivatives(Isigma,'y');
end


%%
% |GRDSMOOTH_CONV|(OLUTION) - The directional gradients' estimation using a
% fast 2D Gaussian convolution of the image with |CONVOLUTION| prior to its
% differentiation with |GRADIENT|.
%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_conv(I, sigma, hsize)                     %#ok

lambda = 3;
if sigma>0.05
    if isempty(hsize)
        hsize = max(1 + round(sigma)*2,7);
    end
    
    % smooth the image out
    [X,Y] = size(I);
    h = gausskernel([hsize hsize], sigma/(lambda*sqrt(X*Y)), [X Y]);
    Ismooth = convolution_base(I, h, 'sym');
else
    Ismooth = I;
end

% calculate the gradient
[gy,gx] = gradient(Ismooth); 
end


%%
% |GRDSMOOTH_VISTA| - The directional gradients' estimation performed in both
% |EDGE| and |CANNYEDGES| functions, using Gaussian filter separability for
% performing directional 1D convolutions.
%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_vista(I, sigma, hsize)                   %#ok

if sigma>0.05
    ssq = sigma^2;
    if isempty(hsize)
        GaussianDieOff = .0001;
        pw = 1:30; % possible hsize
        hsize = find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff,1,'last');
        if isempty(hsize) % still..
            hsize = 1;  % a really small sigma was provided
        end
    end
    
    % design the filters - a Gaussian and its derivative
    t = (-hsize:hsize);
    gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq);     % the gaussian 1D filter
    
    % Find the directional derivative of 2D Gaussian (along I-axis)
    % Since the result is symmetric along I, we can get the derivative along
    % J-axis simply by transposing the result for I direction.
    [x,y] = meshgrid(-hsize:hsize,-hsize:hsize);
    dgau2D = -x.*exp(-(x.*x+y.*y)/(2*ssq))/(pi*ssq);
    
    % Convolve the filters with the image in each direction
    % The canny edge detector first requires convolution with
    % 2D gaussian, and then with the derivitave of a gaussian.
    % Since gaussian filter is separable, for smoothing, we can use
    % two 1D convolutions in order to achieve the effect of convolving
    % with 2D Gaussian.  We convolve along rows and then columns.
    
    %smooth the image out
    ISmooth = imfilter(I,gau,'conv','replicate'); % run the filter across rows
    ISmooth = imfilter(ISmooth,gau','conv','replicate'); % then across columns
    
else
    ISmooth = I;
end

%apply directional derivatives
gy = imfilter(ISmooth, dgau2D, 'conv','replicate');
gx = imfilter(ISmooth, dgau2D', 'conv','replicate');

% Note: using the associativity of convolution and the fact the both the
% Gaussian and the derivative of a Gaussian operations are linear, the first
% derivative of the image function convolved with a Gaussian function is
% equivalent to the image function convolved with the first derivative of
% a Gaussian function. In addition, separability can be used to perform
% convolution on the I and J axes separately.
% Each of the two two-dimensional directional derivatives has to be convolved
% with the image. Using the separability property, two 1D convolutions 
% instead of one costly 2D convolution. That is, the first 1D derivative of
% a Gaussian function convolved with a 1D Gaussian blur function convolved
% with the image function, is equivalent to the image function convolved
% with the first derivative of a 2D Gaussian function. 
end


%%
% |GRDSMOOTH_FLECK| - The directional gradients' estimation implemented
% following the recommandations of [Fleck92] M.M. Fleck: "Some defects in 
% finite-difference edge finders".
%--------------------------------------------------------------------------
function [gv,gh] = grdsmooth_fleck(I, sigma, hsize)                    %#ok               

if sigma>0.05
    if isempty(hsize)
        hsize = round(6*sigma)+1;  % the filter size.
    end
    
    % smooth the image
    gaussian = fspecial('gaussian',hsize,sigma);
    G = filter2(gaussian,I); % imfilter(I, gaussian);
else
    G = I;
end

% differentiate
[rows, cols] = size(I);
h =  [  G(:,2:cols)  zeros(rows,1) ] - [  zeros(rows,1)  G(:,1:cols-1)  ];
v =  [  G(2:rows,:); zeros(1,cols) ] - [  zeros(1,cols); G(1:rows-1,:)  ];
d1 = [  G(2:rows,2:cols) zeros(rows-1,1); zeros(1,cols) ] - ...
    [ zeros(1,cols); zeros(rows-1,1) G(1:rows-1,1:cols-1)  ];
d2 = [  zeros(1,cols); G(1:rows-1,2:cols) zeros(rows-1,1);  ] - ...
    [ zeros(rows-1,1) G(2:rows,1:cols-1); zeros(1,cols)   ];

gh = h + (d1 + d2)/2.0;
gv = v + (d1 - d2)/2.0;
end


%%
% |GRDSMOOTH_LUE(NGO)| - Gaussian derivation recommanded by Luengo using the
% analytical forms of the Gaussian filters and the separability property to
% avoid redundant computations:
%     http://www.cb.uu.se/~cris/blog/index.php/archives/150#more-150
%     http://www.cb.uu.se/~cris/blog/index.php/archives/22
%
% Compare with the method of 'edge': "that code that nicely smoothed the 
% image with a 1D Gaussian first along rows and then along columns, and then
% computed the x and y derivatives by applying full 2D convolutions with 2D 
% Gaussian derivatives [is] wasting a lot of computer clock cycles. Skip the 
% first blurring, [it is not neededd], and compute the gradient as the 
% separable convolution that it is."
%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_luengo(I, sigma, hsize)                   %#ok
[gx,gy] = grdsmooth_lue(I, sigma, hsize);
end % here for convenience (in some calls of the function GRDSMOOTH_BASE)
%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_lue(I, sigma, hsize)                   

if isempty(hsize)
    hsize = 2*ceil(3*sigma)+1;
end
cutoff = floor(hsize/2);

% compute the exact derivative of the Gaussian, and convolve with that
h = fspecial('gaussian',[1,hsize],sigma);
dh = h .* (-cutoff:cutoff) / (-sigma^2);

% I- direction edge detection
gx = conv2(dh,h,I,'same');
% y- direction edge detection
gy = conv2(dh,h,I','same'); gy=gy';
end


%%
% |GRDSMOOTH_ANA|(LYTIC) - The directional gradients' estimation based on an
% analytical form for the Gaussian filter, its derivative and the final 2D
% detector used to convolve the image.
%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_ana(I, sigma, hsize)                      %#ok

if isempty(hsize)
    hsize = round(6*sigma)+1;
end

% J- direction edge detection
fy = d2dgauss(hsize,sigma,hsize,sigma,pi/2);
gy = conv2(I,fy,'same');    % imfilter(I,fy,'replicate','conv');

% I- direction edge detection
fx = d2dgauss(hsize,sigma,hsize,sigma,0);
gx = conv2(I,fx,'same'); 
end


%%
% |D2DGAUSS| - Return a 2D edge detector (first order derivative of a 2D 
% Gaussian function) with size n1*n2; theta is the angle that the detector 
% rotated counter clockwise; and sigma1 and sigma2 are the standard
% deviations of the Gaussian function.
%--------------------------------------------------------------------------
function h = d2dgauss(n1,sigma1,n2,sigma2,theta)
r=[cos(theta) -sin(theta);
   sin(theta)  cos(theta)];

% function for 1D Gaussian filter 
gauss = @(x,std) exp(-x^2/(2*std^2)) / (std*sqrt(2*pi));
dgauss = @(x,std) -x * gauss(x,std) / std^2;                           

h = zeros(n2,n1);
for i = 1 : n2 
    for j = 1 : n1
        u = r * [j-(n1+1)/2 i-(n2+1)/2]';
         h(i,j) = gauss(u(1),sigma1) * dgauss(u(2),sigma2);
        % -u(2) * gauss(u(2),sigma2) * gauss(u(1),sigma1) / sigma2^2;
    end
end

h = h / sqrt(sum(sum(abs(h).*abs(h))));
end


%%
% |GRDSMOOTH_MASK| - The directional gradients' estimation using a fast
% 2D Gaussian convolution of the image with |IMGAUSSIAN| prior to its
% differentiation using directional masks with |GRDMASK|.
%
% See other functions below.
%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_mask(I, sigma, hsize, der)   

if sigma>0.05
    if isempty(hsize)
        hsize = round(6*sigma)+1;   % the filter size.
    end
    
    % smooth the image out
    Isigma = imgaussian(I,sigma,hsize);
else
    Isigma = I;
end

[gx,gy] = grdmask_base(Isigma, der, 'ij');  
end


%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_sobel(I, sigma, hsize)                    %#ok
[gx,gy] = grdsmooth_sob(I, sigma, hsize);
end % here for convenience (in some calls of the function GRDSMOOTH_BASE)
%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_sob(I, sigma, hsize)                      
[gx,gy] = grdsmooth_mask(I, sigma, hsize, 'sobel');
end


%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_prewitt(I, sigma, hsize)                  %#ok
[gx,gy] = grdsmooth_prew(I, sigma, hsize);
end
%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_prew(I, sigma, hsize)                 
[gx,gy] = grdsmooth_mask(I, sigma, hsize, 'prewitt');
end


%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_opt(I, sigma, hsize)                      %#ok
[gx,gy] = grdsmooth_mask(I, sigma, hsize, 'optimal');
end


%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_circ(I, sigma, hsize)                     %#ok
[gx,gy] = grdsmooth_mask(I, sigma, hsize, 'circular');
end

%%
% |GRDSMOOTH_DERIVATIVE5| - The 1st derivatives of the image are estimated 
% using the 5-tap coefficients given by Farid and Simoncelli.  The results 
% are significantly more accurate than MATLAB's |GRADIENT| function on edges
% that are at angles other than vertical or horizontal. 
%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_kovesi(I, sigma, hsize)                   %#ok
[gx,gy] = grdsmooth_tap5(I, sigma, hsize);
end   % for convenience with some calls with 'kovesi' option
%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_derivative5(I, sigma, hsize)              %#ok
[gx,gy] = grdsmooth_tap5(I, sigma, hsize);
end
%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_tap5(I, sigma, hsize)                   
[gx,gy] = grdsmooth_mask(I, sigma, hsize, 'derivative5');
end


%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_derivative7(I, sigma, hsize)              %#ok
[gx,gy] = grdsmooth_tap7(I, sigma, hsize);
end
%--------------------------------------------------------------------------
function [gx,gy] = grdsmooth_tap7(I, sigma, hsize)                 
[gx,gy] = grdsmooth_mask(I, sigma, hsize, 'derivative7');
end


