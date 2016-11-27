%% FASTCROSSDILATION - NOT FINISHED IMPLEMENTING
%
%% Algorithm
% # inputs:
%
%   I : image I(r) 
%   z : random noise z(r) 
%   Vx, Vy : x and y components of v(r) 
%   niter : number of iterations 
%   h : step size of the Euler algorithm 
%
% # output:
%
%   I1 : filtered image
%
% # pseudo-code
%
%   for i = 1 : N iter  { 
%      I1 = evaluate(I, Vx, Vy, h); 
%      z1 = evaluate(z, Vx, Vy, h); 
%      for each (x,y)  { 
%         if z(x,y) ? z1(x,y) 
%         I(x,y) = I1(x,y); 
%         } 
%      z = max(z, z1); 
%      Wx = evaluate(Vx, Vx, Vy, h); 
%      Wy = evaluate(Vy, Vx, Vy, h); 
%      Vx = Vx + Wx; 
%      Vy = Vy + Wy; 
%   } 
%
%% Reference
% [PP09]  G. Papari and N. Petkov: "Spatially variant dilation for 
%      unsupervised painterly rendering", Proc. ISMM, 2009. 

%% Function implementation
function U = fastcrossdilation(I, z, v, niter, h) 

t = isequal(z,I);

for i=1:niter
    U = evaluate(I, v(:,:,1), v(:,:,2), h);
    if t,  z1 = U;
    else   z1 = evaluate(z, v(:,:,1), v(:,:,2), h);   end
    I(z<=z1) = U(z<=z1);
    z = max(cat(3,z,z1),[],3);
    wx = evaluate(v(:,:,1), v(:,:,1), v(:,:,2), h); 
    wy = evaluate(v(:,:,2), v(:,:,1), v(:,:,2), h); 
    v(:,:,1) = v(:,:,1) + wx;
    v(:,:,2) = v(:,:,2) + wy;
end

end % end of fastcrossdilation


%%
% |EVALUATE| - Take in input three real valued functions $I(x,y)$, $V_x(x,y)$,
% and $V_y(x,y)$, with $(x,y) \in Z^2$, and a scalar $h$, and returns the
% function $U(x,y) = I[x+ h \cdot V_x(x, y), y + h \cdot V_y(x,y)]$, which is
% computed by means of bilinear interpolation
% -------------------------------------------------------------------------
function a1 = evaluate(a, vx, vy, h)

[X,Y] = size(a);
[x,y] = meshgrid(1:X,1:Y);
xi = x + h*vx;
yi = y + h*vy;
a1 = interp2(a, xi, yi, 'linear');

end % end of evaluate

