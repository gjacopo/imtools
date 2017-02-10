function rec = singfractalrecons_base(msm, gx, gy)
% SINGFRACTALRECONS_BASE.
%
%       rec = SINGFRACTALRECONS_BASE(msm, gx, gy);
%
% credit: J.Grazzini (ISR-2/LANL)
%
% See also 
% Related: SINGFRACTAL_BASE, FRACTALWAVE_BASE.

rec = reconstruct_naif (gx(msm), gy(msm));
end


%--------------------------------------------------------------------------
function rec = reconstruct_naif (gx, gy)

[X Y] = size(gx);

% Direct Fourier Transform of the signal in Gx
FFTgx = fft2(gx);
% Direct Fourier Transform of the signal in Gy
FFTgy = fft2(gy);

[x y] = meshgrid(0:X-1, 0:Y-1);
I = x>X/2; x(I) = x(I) - X;
I = y>Y/2; y(I) = y(I) - Y;
% [X Y] = meshgrid( [0:Y/2-1,0,-Y/2+1:-1]/Y, ...
% 		  [0:X/2-1,0,-X/2+1:-1]/X);
% frequency vector
dX = sin(pi.*x);
dY = sin(pi.*y);
% or: dX =X; dY =Y;

prefm = dX.^2 + dY.^2;
prefm = (prefm>1e-30) .* prefm + (prefm<=1e-30).*1;

% The scalar product with the frequency vector f is computed and the result
% is multiplied by the imaginary unit i=sqrt(-1), wich correspond to a 
% rotation:
%     A = a_0 + a_1*i   =>   B = A*i = -a_1 +a_0*i         
gxR = (dX.*imag(FFTgx) + dY.*imag(FFTgy)) ./ prefm;
gxI = - (dX.*real(FFTgx) + dY.*real(FFTgy)) ./ prefm;

% Inverse Fourier transform 
rec = real(ifft2( gxR + 1i*gxI ));

end