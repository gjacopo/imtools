%% REDUCE2D_BASE - Base function for REDUCE2D.
%
%% Syntax
%    Ir = REDUCE2D_BASE(I, filter, order);
%
%% See also
% Related:
% <REDUCE2D.html |REDUCE2D|>, 
% <PYRDOWN.html |PYRDOWN|>, 
% <EXPAND2D_BASE.html |EXPAND2D_BASE|>.
% Called: 
% |REDUCE2D_MEX|.

%% Function implementation
function Ir = reduce2d_base(I, filter, order)

%% 
% dealing with multispectral images
[X,Y,C] = size(I);
if X<2 || Y<2,  return;  end;   % otherwise we can crash matlab!!!

if C>1
    Ir = zeros(floor(X/2),floor(Y/2),C);
    for c=1:C
            Ir(:,:,c) = reduce2d_base(I(:,:,c), filter, order);
    end
    return;
end

%%
% main program

if strcmpi(filter,'lpl')
    %%
    % classical Laplacian pyramid: Matlab implementation
    h = [1 4 6 4 1]/16;
    Ir = conv2(h, h, I, 'same');
    Ir = Ir(1:2:X,1:2:Y);

else
    %%
    % calling the mex program for spline based pyramids
    Ir = reduce2d_mex(I, filter, order);
end

end % end of reduce2d_base