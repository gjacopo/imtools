%% EXPAND2D_BASE - Base function for EXPAND2D. 
%
%% Syntax
%    Ie = EXPAND2D_BASE(I, filter, order);
%
%% See also
% Related:
% <EXPAND2D.html |EXPAND2D|>, 
% <PYRUP.html |PYRUP|>,
% <REDUCE2D_BASE.html |REDUCE2D_BASE|>. 
% Called: 
% |EXPAND2D_MEX|.
                                                    
%% Function implementation
function Ie = expand2d_base(I, filter, order)

%% 
% dealing with multispectral images
[X,Y,C] = size(I);

if C>1
    Ie = zeros(X*2,Y*2,C);
    for c=1:C
        Ie(:,:,c) = expand2d_base(I(:,:,c), filter, order);
    end
    return;
end

%%
% main program

if strcmpi(filter,'lpl')
    %%
    % classical Laplacian pyramid: Matlab implementation
    h = [1 4 6 4 1]/16;
    Ie = zeros(2*X, 2*Y);
    Ie(1:2:2*X,1:2:2*Y) = I;
    Ie = 4 * conv2(h, h, Ie, 'same');
    
else
    %%
    % calling the mex program for spline based pyramids
    Ie = expand2d_mex(I, filter, order);
end

end % end of expand2d_base
