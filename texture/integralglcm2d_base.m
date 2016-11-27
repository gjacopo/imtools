%% INTEGRALGLCM2D_BASE - Base function for INTEGRALGLCM2D.
%
%% Syntax
%     O = INTEGRALGLCM2D_BASE(I, res, dcar, win);
%
%% See also
% Related:
% <INTEGRALGLOV2D_BASE.html |INTEGRALGLOV2D_BASE|>,
% <LOCALGLCM2D_BASE.html |LOCALGLCM2D_BASE|>,
% <INTEGRALIMAGE.html |INTEGRALIMAGE|>.
% Called:
% <SLIDEHISTOFUN_BASE.html |SLIDEHISTOFUN_BASE|>.

%% Function implementation
function O = integralglcm2d_base(I, res, dcar, win )

%%
% checking parameters and setting internal variables

[X,Y] = size(I(:,:,1));                                                    

imax = max(I(:)); imin = min(I(:));
if nargin<4,   win = 3;
    if nargin<3,  dcar = [1 1];
        if nargin<2,  res = round(imax - imin) + 1;
        end
    end
end

if imax==imin,
    error('integralglcm2d_base:inputerror', 'constant image - nothing to do');
elseif res==0
    res = round(imax - imin) + 1; % again, default
end

%%
% quantize the input image
I = floor((res-1) * (I - imin) / (imax - imin));

%%
% ensure that the window size is odd
if mod(win,2)==0,  win = win + 1;  end

%%
% use the Cartesian representation in the following
dnorm = sqrt(sum(dcar.^2,2));         
nd = length(dnorm); % also: size(dcar,1);

% useful indexes
% pad = ceil(max(dnorm));
% pixI = reshape(1:X*Y,[X Y]);
% index of the displaced images
% pixT = padarray(pixI, [pad pad], 'replicate', 'both');
% pixIinT = reshape(pad*((X+1):(X+Y)), [X 1]) + pixI;

O = cell(nd,1);

%%
% handle of the texture feature functions
func  = {@contrast, @energy, @homogeneity, @entropy2};

%% 
% main computation

for d=1:nd
    dx = dcar(d,1); dy = dcar(d,2);
    Id = padarray(I, [abs(dx) abs(dy)], 'symmetric', 'both' );
    Id = Id(abs(dx)+dx+1:abs(dx)+dx+X,abs(dy)+dy+1:abs(dy)+dy+Y,:);
%    O{d} = slidehistofun_base(I + res * Id, 'dist', func, win, res^2);
    O{d} = slidehistofun_base(I + res * Id, func, win, res^2, 'dist');
    % O = slidehistofun_base(I + res * Id, 'int', @contrast, win, res^2);
end

if nd==1,  O = O{1};  end

end % end of integralglcm2d_base


%% Subfunctions

%--------------------------------------------------------------------------
function E = energy(pI)                                               
E = sum(pI.^2,1);
end % end of energy


%--------------------------------------------------------------------------
function E = entropy2(pI)                                              
pI(pI==0) = 1;
E = - sum(pI .* log2(pI),1);
end % end of entropy2


%--------------------------------------------------------------------------
function H = homogeneity(pI)                                           
res = sqrt(size(pI,1));
[I,J] = meshgrid(1:res,1:res);
H = 1 + (I-J).^2;
H = repmat(H(:),[1 size(pI,2)]);
H = sum(pI ./ H, 1);
end % end of homogeneity


%--------------------------------------------------------------------------
function C = contrast(pI)                                             
res = sqrt(size(pI,1));
[I,J] = meshgrid(1:res,1:res);
C = (I-J).^2;
C = repmat(C(:),[1 size(pI,2)]); 
C = sum(pI .* C,1);
end % end of contrast
