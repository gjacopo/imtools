%% INTEGRALGLOV2D_BASE - Base function for INTEGRALGLOV2D.
%
%% Syntax
%     O = INTEGRALGLOV2D_BASE(I, res, win);
%
%% See also
% Related:
% <INTEGRALGLCM2D_BASE.html |INTEGRALGLCM2D_BASE|>,
% <LOCALGLOV2D_BASE.html |LOCALGLOV2D_BASE|>,
% <INTEGRALIMAGE.html |INTEGRALIMAGE|>.
% Called:
% <SLIDEHISTOFUN_BASE.html |SLIDEHISTOFUN_BASE|>.

%% Function implementation
function O = integralglov2d_base(I, res, win)

%%
% checking parameters and setting internal variables

imax = max(I(:)); imin = min(I(:));
if nargin<3,   win = 3;
        if nargin<2,  res = round(imax - imin) + 1;
    end
end

if imax==imin,  
    error('integralglov2d_base:inputerror', 'constant image - nothing to do');
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
% handle of the texture feature functions
func  = {@contrast, @energy, @homogeneity, @entropy2, @variance};

%%
% main computation

%O = slidehistofun_base(I, 'int', @entropy2, win, res);
O = slidehistofun_base(I, 'dist', func, win, res);
%O = slidehistofun_base(I, 'dist', @energy, win, res);

end % end of integralglov2d_base


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
H = 1 + repmat((1:size(pI,1))', [1 size(pI,2)]).^2;
H = sum(pI ./ H, 1); 
end % end of homogeneity


%--------------------------------------------------------------------------
function C = contrast(pI)                                         
C = repmat((1:size(pI,1))', [1 size(pI,2)]).^2;
C = sum(pI .* C,1); 
end % end of contrast


%--------------------------------------------------------------------------
function V = variance(pI)                                         
I = repmat((1:size(pI,1))', [1 size(pI,2)]);
muI = repmat(sum(I .* pI,1), [size(pI,1) 1]);

V = sum (pI .* (I .^2) - muI.^2,1);
end % end of variance
