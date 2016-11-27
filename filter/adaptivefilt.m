%% ADAPTIVEFILT - Perform adaptive filtering.
%
%% Syntax
%     F = ADAPTIVEFILT(I, H, M);
%
%% Inputs
% *|I|* : a |(n1,n2)| input image.
%
% *|H|* : array of |m| kernels: it is a |(p1,p2,m)| matrix where each |H(:,:,k)| 
%     is a filter; |p1| and |p2| should be odd integers; note that during 
%     filtering, |H(:,:,k)| is automatically normalized to 1.
%
% *|m|* : a |(n1,n2)| matrix of integer in |{1,...,m}|, where |m(i,j)| is the 
%     index of the kernel to be used in pixel |(i,j)|.
%
%% Output
% *|F|* : a |(n1,n2)| filtered image.
%
%% See also
% Related:
% <CONVOLUTION.html |CONVOLUTION|>,
% <TENSANIFILT.html |TENSANIFILT|>,
% <TENSCALEDIRFILT.html |TENSCALEDIRFILT|>,
% <GEODESICFILT.html |GEODESICFILT|>,
% <MDLFILT.html |MDLFILT|>.
% Called:
% <ADAPTIVEFILT_BASE.html |ADAPTIVEFILT_BASE|>.

%% Function implementation
function F = adaptivefilt(I, H, m)

%% 
% check if possible
if ~exist('adaptivefilt_mex','file')
    error('adaptivefilt:missinglibrary', ...
        'mex file adaptivefilt_mex not found');
end

error(nargchk(1, 3, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

if ~isnumeric(I)
    error('adaptivefilt:inputparameter','a matrix is required in input'); 
end

%%
% main calculation

F = adaptivefilt_base(I, H, m);

%%
% display

if p.disp
    figure, imagesc(rescale(F,0,1)), axis image off;
    title('Adaptively filtered image');  if size(F,3)==1, colormap gray;  end;
end

end % end of adaptivefilt
