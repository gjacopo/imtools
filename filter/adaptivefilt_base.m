%% ADAPTIVEFILT_BASE - Base function for ADAPTIVEFILT.
%
%% Syntax
%   F = ADAPTIVEFILT_BASE(A, H, I);
%
%% See also
% Related:
% <ADAPTIVEFILT.html |ADAPTIVEFILT|>,
% <TENSCALEDIRFILT_BASE.html |TENSCALEDIRFILT_BASE|>,
% <TENSANIFILT_BASE.html |TENSANIFILT_BASE|>,
% <CONVOLUTION_BASE.html |CONVOLUTION_BASE|>,
% <GEODESICFILT_BASE.html |GEODESICFILT_BASE|>,
% <MDLFILT_BASE.html |MDLFILT_BASE|>.
% Called:
% ADAPTIVEFILT_MEX.

%% Function implementation
function F = adaptivefilt_base(I, H, M)

%% 
% dealing with multispectral images
C = size(I,3);

if C>1
    F = I;
    for c=1:C
        F(:,:,c) = adaptivefilt_base(I(:,:,c), H, M);
    end
    return;
end

%% 
% call the mex file

F = adaptivefilt_mex(I, H, M);

end % end of adaptivefilt_base
