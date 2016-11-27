function [ZCmap,GradSq] = sarkaredge(I, alpha, beta,Imult,ZCThr)
% SARKAREDGE - Implements the Optimal Zero Crossing edge Operator (OZCO)

% gamma = scale;
% alpha in -1.5 to 1.5, the OZCO shape parameter; default alpha=0.312
% beta in 0.2 to 2.0, the OZCO scale parameter; default: beta=1/gamma
	% Imult : image brightness multiplier for ZC in ROI 
	% ZCThr : threshold for ZC validation 
%   [SB91] 	S. Sarkar and K. L. Boyer: "Optimal Infinite Impulse response zero
%      crossing based edge detectors", CVGIP, 54(2):224-243, 1992.


%I = I0/max(I0(:));
[sizeR,sizeC] = size(I);

if exist('sarkaredge_mex','file')
    warning('sarkaredge_base:error', ...
        ['TODO:wrong version of sarkaredge_mex implemented - '...
        'currently using matlab code']);
end   


%% Build the causal and non causal filters
% 	- px (py) is the smoothing function (second integral of OZCO).
% 	- dx (dy) is the smoothed differentiator (first integral of OZCO).

% the filter is performed in 2D by applying the optimal filter perpendicular
% H to the edge and a projection function G parallel to the edge
% i) apply the noncausal projection filter G(z) to the rows of the im0.3age and
%  the noncausal edge detection filter H(z) to the columns; this estimates
%  the smoothed image gradient in the row direction. 
% ii) similarly, apply the noncausal version of G(z) to the columns and the
%  noncausal version of H(z) to the rows to get an estimate of the smoothed
%  gradient in the column direction. 
% The number of multiplications per pixel is 40 and is independent of the 
% size (scale) of the filter used.

% 	First, some constants needed...
% sampling interval
tau = 1; % tau = [taux tauy];

alpha0 = (alpha+2) / beta^2;
alpha1 = exp(-beta*tau) * ...
    (tau*(alpha+2)/beta + 0.5*tau^2*(alpha+1) - 2*(alpha+2)/beta^2);
alpha2 = exp(-2*beta*tau) * ...
    (-tau*(alpha+2)/beta + 0.5*tau^2*(alpha+1) + (alpha+2)/beta^2);

beta1 = -3 * exp(-beta*tau);
beta2 = 3 * exp(-2*beta*tau);
beta3 = - exp(-3*beta*tau);

gamma1 = (0.5*beta*tau^2*(alpha+1) + tau) * exp(-beta*tau);  %a01
gamma2 = (0.5*beta*tau^2*(alpha+1) - tau) * exp(-2*beta*tau); % a02

nconst = 2 * (3*alpha+5) / beta^3;

% expression of the projection IIR filter:
%     G(z) = (a0 + a1*z + a2*z^2) / (1 - b1 *z - b2*z^2 - b3*z^3)
% the projection function is chosen to be the integral of the edge detection
% filter, which is a low-pass filter very similar to a Gaussian
% this noncausal filter can be realized by the sum of two identical causal
% filters operating in opposite directions
% see functions derivy_convl and derivx_convl

% build the smoothing function feedforward (causal) coefficient vector
pb = [alpha0 alpha1 alpha2];  % used with half_projec_convl
pb = pb * beta^2;

% build the smoothing function feedback (anti causal) coefficient vector
pa = [1 beta1 beta2 beta3];

% build the differentiator feedforward coefficient vector
db = [0 gamma1 gamma2]; % used with half_opt_convl

% the differentiator feedback coefficients are the same as the smoother's
da = pa; 

% OZCO Infinite Impulse Responses: run the filters over a unit pulse and
% display the response
%test = [zeros(1,100) 1 zeros(1,100)];
%pxr = filter(pb,pa,test);
%dxr = filter(db,da,test);
%pxl(201:-1:1) = filter(pb,pa,test(201:-1:1));
%dxl(201:-1:1) = filter(db,da,test(201:-1:1));
%px = pxl + pxr - [zeros(1,100) pxr(101) zeros(1,100)];
%dx = dxr - dxl;
%plot(px,'g')
%plot(10*dx,'b')

% Obtain zc detection and validation parameters 
ZCThr = -ZCThr;
%	ImThr = Imult*max(max(I));

%% Apply filters
% convolve with differentiator in X, smoothing function in Y for the X-
% component of the smoothed gradient.
% first apply opt_convl
Ifx = filter(db, da, I, [], 2);
Ibx = fliplr(filter(db, da, fliplr(I), [], 2)); 
Ifx = Ifx - Ibx;   
% then apply projec_convl 
Ibx = flipud(filter(pb, pa, flipud(Ifx), [], 2));
Ifx = filter(pb, pa, Ifx, [], 1);
Ifx = Ifx + Ibx - alpha0*Ifx;

% convolve with differentiator in Y, smoothing Function in X for the Y-
% component of the smoothed gradient.
% first apply projec_convl
Ify = filter(pb, pa, I, [], 2);
Iby = fliplr(filter(pb, pa, fliplr(I), [], 2));
Ify = Ify + Iby - alpha0*I;
% then apply opt_convl 
Iby = flipud(filter(db, da, flipud(Ify), [], 2));
Ify = filter(db, da, Ify, [], 1);
Ify = Ify - Iby;

%% Inner product of gradient with itself for the gradient squared magnitude.
Ig = (Ifx .* Ifx) + (Ify .* Ify);
max(Ig(:))
Itheta = atan2(Ify,Ifx);

% compute the gradient of the gradient squared magnitude for the first
% directional derivatives
[gx,gy] = gradient(Ig); 

% compute the inner product of the gradient of gradient squared with the
% original smoothed gradient for the second directional derivative in the
% direction	of the smoothed gradient.  
% zero crossings in this image indicate potential edges, subject to
% qualification and strength tests.
DI = Ifx .* gx + Ify .* gy;
DI = DI ./ (2 * tau * Ig);
DI = DI / nconst;
DI(Ig==0) = 0;
max(DI(:))



%% Detect and qualify zcs

% ZCs equal border pixels in DI; corresponding Image pixel > ImThr
ZCmap = zeros(sizeR,sizeC);
% ZC2 = zeros(sizeR,sizeC);
% DxI= zeros(sizeR,sizeC);
% GradX = zeros(sizeR,sizeC);
% GradY = zeros(sizeR,sizeC);

% computing initial ZC map
r = 2:sizeR-1;
c = 2:sizeC-1;
ImThr = Imult*max(Ig(:));

%ZCmap(r,c) = I(r,c)>ImThr & DI(r,c) >= 0 & (DI(r-1,c)<0 | ...
%   DI(r+1,c)<0 |DI(r,c+1)<0 | DI(r,c-1)<0);

size(Ig(r,c))
max(max(Ig(r,c)))
figure, imagesc(Ig(r,c)>ImThr)
ZCmap(r,c) = Ig(r,c)>ImThr & DI(r,c) >= 0 & ...
    (DI(r-1,c)<0 | DI(r+1,c)<0 |DI(r,c+1)<0 | DI(r,c-1)<0);
max(ZCmap(:))

figure, imagesc(ZCmap)
% smoothed gradient consistent with third directional derivative
% compute normalized third directional derivatives
[GradX,GradY] = gradient(DI);
DxI = (Ifx .* GradX) + (Ify .* GradY);
DxI = DxI/max(max(abs(DxI)));

% pruning weak ZCs from map

ZCmap = ZCmap & (DxI < ZCThr);

% pruning ZCs inconsistent with gradient orientation
% DxI = (DxI .* GradSq);
% ZCmap = ZCmap & (DxI < 0);

% qualified ZCmap : phantoms deleted
%    1-ZCmap

% qualified ZCmap overlay : strong ZCs
%    max(ZCmap,Image)
end

