%% ANISORFILT_BASE - Base function for ANISORFILT.
%
%% Syntax
%     u = ANISORFILT_BASE(I, delta, niter, nu, rho, sigma, der, int, samp);
%
%% See also
% Related:
% <ANISORFILT.html |ANISORFILT|>.
% Called:
% <ANISOR_BASE.html |ANISOR_BASE|>.

%% Function implementation
function u = anisorfilt_base(I, niter, delta, nu, rho, sigma, der, int, samp)

%%
% setting internal variables

[X,Y] = size(I(:,:,1)); 

% predefine appropriate indices
zer0 = false(X+2,Y+2);

iminusj = zer0; iminusj(1:end-2,2:end-1) = true; iminusj = iminusj(:);
iplusj = zer0; iplusj(3:end,2:end-1) = true; iplusj = iplusj(:);
ijminus = zer0; ijminus(2:end-1,1:end-2) = true; ijminus = ijminus(:);
ijplus = zer0; ijplus(2:end-1,3:end) = true; ijplus = ijplus(:);

zer0 = zeros(X*Y,1);

%%
% perform Rouy-Turin numerical scheme like in Eq.(11)

u = double(I);

for i =1:niter
    kappa = anisor_base(u, nu, rho, sigma, der, int, samp);
    %  kappa=1;

    U = padarray( u, [1 1], 'replicate', 'both');
    
    % Rouy-Tourin scheme a in Eq.(11): naive implementation
    h = 1;
    u = u(:) + delta * kappa(:) .* ...
        sqrt((max([zer0 (U(iplusj) - u(:))/h (U(iminusj) - u(:))/h],[],2)).^2 + ...
        (max([zer0 (U(ijplus) - u(:))/h (U(ijminus) - u(:))/h],[],2)).^2);
    % note: the sign of delta defines if the morphological operation is a
    % dilation ('+' sign) or an erosion ('-' sign)
    
    % http://cermics.enpc.fr/~forcadel/Publi/FGL.pdf
    % http://etd.lsu.edu/docs/available/etd-09152008-143521/unrestricted/Tu
    % gurlandiss.pdf
    % http://ctr.stanford.edu/ResBriefs03/herrmann1.pdf
    % http://www.mia.uni-saarland.de/Publications/pizarro-ismm09.pdf
   
    u = reshape(u,[X,Y]);
end

end % end of anisorfilt_base
