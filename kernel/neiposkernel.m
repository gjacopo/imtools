%% NEIPOSKERNEL - Inscribed window indexing.
%
%% Description
% Create the kernel containing the index (and also subscript) positions of
% the grid pixels embedded within a square window of given radius within an 
% image with given number of rows (I-dimension) relative to the center of
% this window. 
%
%% Syntax
%     ind = NEIPOSKERNEL(ws, X);
%     [ind, ix, iy] = NEIPOSKERNEL(ws, X);
% 
%% Inputs
% *|ws|* : radius of the square window.
%
% *|X|* : number of rows. 
% 
%% Outputs 
% *|ind|* : index positions relative to the central pixel of the window of 
%     radius ws of the pixels contained in this window.
%
% *|ix, iy|* : subscript position relative to the central pixel (nothing else
%     than the outputs of |NDGRID|).
%
%% See also
% Called:
% <matlab:webpub(whichpath('NDGRID')) |NDGRID|>.

%% Function implementation
function [ind, varargout] = neiposkernel(ws, X)

[ix, iy] = ndgrid(-ws:ws);
ind = ix + X*iy;
% ind = ind(:);

if nargout>1,  varargout{1} = ix;
    if nargout>2,  varargout{2} = iy;  end
end

end % end of neiposkernel
