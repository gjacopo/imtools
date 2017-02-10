function Tpsi = waveprofile_base(I, rho, measure, sigma, der, int, samp, gn, tn)
% WAVEPROFILE_BASE - Perform a wavelet decomposition.
%
%    Tpsi = WAVEPROFILE_BASE(I, rho, measure, sigma, der, int, samp, gn, tn);
%
% credit: J.Grazzini (ISR-2/LANL)
% 
% See also 
% Related: FRACTALWAVE_BASE --
% Called: GSTSMOOTH_BASE, GSTDECOMP, IFFT2, FFT2.

%funeig = @(x,y)sqrt(x+y);  % trace: edges
%funeig = @(x,y)((x.*y)./(x+y));  % det/trace: corners
%funeig = @(x,y)x;
funeig = @(x,y)sqrt(x.*y);  % sqrt(det)

n_rho =  length(rho);

% array storing the wavelet transforms at different scales
Tpsi = [];

if any(strcmp(measure,{'grad','scalar','decorr'}))
    % compute the gradient norm once for all
    %[gx,gy] = gradient(I);
    %mag = sqrt(gx.^2+gy.^2);
    mag = modgradient( I );
end

for r=1:n_rho
    
    % current scale
    rr = rho(r);
    
    if any(strcmp(measure,{'grad','scalar','decorr'}))
        
        % % compute the wavelet transform using a Gaussian function T_psi
        % % note that computation are made in the spatial domain
        % wave = fspecial('gaussian',[1,h],rho);
        % % compute the wavelet transform by line- then column 1D convolutions
        % Tpsi_r = conv2(wave,wave,mag,'same') / sqrt(XY);
        wave = wavegauss(rr, [size(I,1), size(I,2)]);
        % Tpsi_r = ifftshift(real(ifft2(fft2(fftshift(wave)).*fft2(mag))));
        Tpsi_r = real(ifft2(fft2(wave).*fft2(mag)));
        
    elseif any(strcmp(measure,{'tens','riemann','corr'}))
        % compute the GST at given scales
        T = gstsmooth_base(I, rr, sigma, der, int, samp, [], gn, tn, [], []);
        % T = T(h+1:X+h, h+1:h+Y,:,:);
        
        % estimate the spatial average \bar{T}(rho)
        %meanT = mean(reshape(T,[n*m 4]),1);
        %meanT = reshape(ones(n*m,1)*meanT,[n m 2 2]);
        % transform :
        %   T_{\Psi} s(\vec{x},rho) / \bar{T}(rho))
        %T = T ./ meanT;
        % normalize by the dimensions of the image
        T = T/numel(T(:,:,1,1));
        
        [Tpsi_r,Tpsi_r2] = gstdecomp(T);
        
        Tpsi_r = funeig(Tpsi_r,Tpsi_r2); % first eigenvalue
    end
    
    % concatenate the current results with those computed at the previous
    % scales
    Tpsi = cat(3, Tpsi, Tpsi_r);
    % Tpsi is a 3D array (X x Y x nsigma), where the transform at
    % scale qs^i are stored in the i-th frame
    
end

end


%--------------------------------------------------------------------------
function [psi,sc] = wavegauss( scale, N )

X = N(1); Y = N(2);
% [X Y] = meshgrid( [0:X/2-1,-X/2:-1], ...
% 		  [0:Y/2-1,-Y/2:-1] )
[x y] = meshgrid(0:X-1, 0:Y-1);
I = x>X/2; x(I) = x(I) - X;
I = y>Y/2; y(I) = y(I) - Y;

%[X Y] = meshgrid( [-X/2:-1, 0:X/2-1], ...
%    [-Y/2:-1,0:Y/2-1]  );
% the previous meshgrid avoids the futher call to fftshift
% in convolution processes...

R = (x.*x + y.*y);
psi = exp(- R / (2*scale^2)); % / scale^2;

s=sum(abs(psi(:)));
if s>1.e-17
    psi = psi/s;
end
norma = s / (X*Y);
sc = sqrt(norma/pi);
end


%--------------------------------------------------------------------------
function [mag, ax, ay] = modgradient( f )

j=sqrt(-1);

[X Y] = size( f );

g = f;
% direct Fourier transform of the signal in gx
FFTgx = fft2(g);
% which is equivalent to: FFTgx = fft(fft(gx/sqrt(xeff)).'/sqrt(yeff)).';

% [X Y] = meshgrid( [0:yeff/2-1,0,-yeff/2+1:-1]/yeff, ...
% 		  [0:xeff/2-1,0,-xeff/2+1:-1]/xeff);
[x y] = meshgrid(0:X-1, 0:Y-1);
I = x>X/2; x(I) = x(I) - X;  x = x / X;
I = y>Y/2; y(I) = y(I) - Y;  y = y / Y;

% frequency vector
dx = 2 * sin(pi.*x);
dy = 2 * sin(pi.*y);
% or: dX =X; dY =Y;

% inverse Fourier transform
ax = real(ifft2( - dx.*imag(FFTgx) +  j*dx.*real(FFTgx) ));
ay = real(ifft2( - dy.*imag(FFTgx) +  j*dy.*real(FFTgx) ));

mag = sqrt(ax.*ax + ay.*ay);

end


