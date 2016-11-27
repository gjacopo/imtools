%% GRDMASK - Mask-based directional derivatives of an image.
%
%% Description
% Apply classical gradient operators using local difference masks.
%
%% Syntax
%    [Gx,Gy] = GRDMASK(I);
%    [Gx,Gy] = GRDMASK(I, method, map, axis);
%
%% Inputs
% *|I|* : input 2D image with size |(X,Y,C)|, with |C>1| for multispectral image.
%
% *|method|* : (optional) parameter defining the original |(3,3)| gradient 
%      mask used for estimating the directional derivatives [GW02]; classical
%      masks include:
% * |'sob'| or |'sobel'| for Sobel masks in i- and j-directions (see option
%          |axis| below),
% * |'kir'| or |'kirsch'| for Kirsch mask,
% * |'prew'| or |'prewitt'| for Prewitt mask,
% * |'iso'| or |'isotropic'| for Frei-Chen's isotropic mask,
% * |'opt'| or |'optimal'| for Ando's optimal mask,
% * |'ori'| or |'orientation'| for optimized orientation invariant filter,
% * |'circ'| or |'circular'| for Davies' circular mask,
% * |'rob'| or |'robinson'| for Robinson mask,
%
% finite difference schemes are also implemented as:
%
% * |'matlab', 'diff'| or |'difference'|  when central differences are
%          computed (using Matlab function |GRADIENT|),
% * |'for'| or |'forward'|  for forward differences,
% * |'back'| or |'kackward'|  for backward differences,
% 
% but are not very accurate when further used for the detection of edges 
%      that are at angles other than vertical or horizontal [Fleck92],
%      therefore improvements have been proposed in [FS04] and are here
%      implemented through the call to P.Kovesi's functions [KOVESI] as:
%
% * |'tap5'| or |'derivative5'| using 5-tap coefficients,
% * |'tap7'| or |'derivative7'| using 7-tap coefficients;
%
%      default: |method='sob'|. 
%
% *|map|* : (optional) logical map (with size |(X,Y)|) defining the pixels of
%      the input image domain where the gradient is to be computed; set to
%      false where no gradient needs to be output; in this case, the function
%      |GRDMASKMAP_BASE| is called; default: |map=[]|, ie. considered to be
%      true everywhere and the function |GRDMASK_BASE| is called.
%
% *|axis|* : (optional) string indicating the direction/orientation of the output
%      directional derivatives, it is either |'ij'| or |'xy'|; default: |axis='ij'|
%      (right-handed coordinate system), ie. the vector orthogonal to the
%      (real) gradient is output: |gx| is in that case the vertical derivative,
%      and |gy|, the horizontal derivative.
%
%% Outputs
% *|Gx, Gy|* : directional derivatives in vertical and horizontal directions 
%      resp. in the case |axis='ij'|, the opposite in the case |axis='xy'|.
%
%% References
% [Fleck92] M.M. Fleck: "Some defects in finite-difference edge finders",
%      IEEE Trans. Pattern Analysis and Machine Intelligence, 14(3):337-345,
%      1992.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=120328&tag=1>
%
% [GW02]  R.C. Gonzales and R.E. Woods: "Digital Image Processing", Prentice     
%      Hall, 2002.
%
% [FS04]  H. Farid and E. Simoncelli: "Differentiation of discrete multi-
%      dimensional signals", IEEE Trans. Image Processing, 13(4):496-508,
%      2004.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1284386>
%
% [KOVESI]  P.D. Kovesi: "MATLAB and Octave Functions for Computer Vision
%      and Image Processing", The University of Western Australia, available
%      at <http://www.csse.uwa.edu.au/~pk/research/matlabfns/>.
%
%% Acknowledgments: 
% Options |'derivative5'| and |'derivative7'| call the resp. functions of
% P.Kovesi's toolboox available at:
%        http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/.
%
%% See also  
% Related:
% <GRDSMOOTH.html |GRDSMOOTH|>,
% <matlab:webpub(whichpath('GRADIENT')) |GRADIENT|>.
% Called:
% <GRDMASK_BASE.html |GRDMASK_BASE|>,
% <GRDMASKMAP_BASE.html |GRDMASKMAP_BASE|>.

%% Function implementation
function [Gx,Gy] = grdmask(I, varargin)

%% 
% parsing parameters 

error(nargchk(1, 12, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

if ~isnumeric(I)
    error('grdmask:inputerror','a matrix is required in input'); 
end

p = createParser('GRDMASK');   
p.addOptional('method', 'sobel', @(x)ischar(x) && ...
    any(strcmpi(x,{'matlab','diff','difference', 'backward','back', 'forward','for',...
    'prewitt','prew', 'kirsch','kir', 'robinson','rob', 'circular','circ',...
    'optimal','opt', 'ori','orientation', 'isotropic','iso', 'sobel','sob', ...
    'roberts', 'derivative5','tap5', 'derivative7','tap7'})));
p.addOptional('map',[],@(x)isnumeric(x) || islogical(x));
p.addOptional('axis','ij',@(x)ischar(x) && any(strcmpi(x,{'ij','xy'})));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p); 

%% 
% check parameters

if any(strcmpi(p.method,{'derivative5','tap5', 'derivative7','tap7'}))
    error('grdmask:methoderror', ...
        'Kovesi''s library required with tap coefficients');

elseif ~isempty(p.map) && ~isequal(size(p.map),size(I(:,:,1)))
    error('grdmask:inputerror', ...
        'input map must be same size as the input image domain');
end

%% 
% main computation

[X,Y,C] = size(I);

if isempty(p.map) || all(p.map(:))
    [Gx,Gy] = grdmask_base(I, p.method, p.axis);
else
    [Gx,Gy] = grdmaskmap_base(I, p.map, p.method, p.axis);
end

%%
% display

if p.disp
    figure, 
    subplot(1,2,1), imagesc(rescale(Gx,0,1)), axis image off, 
    title(['gx - axis ''' p.axis '''']); if C==1,  colormap gray;  end
    subplot(1,2,2), imagesc(rescale(Gy,0,1)), axis image off, 
    title(['gy - axis ''' p.axis '''']); if C==1,  colormap gray;  end
    figure,
    imagesc(rescale(I,0,1)), axis image off;    hold on;
    if C==1,  colormap gray;  end
    [X,Y] = meshgrid(3:10:X-2,3:10:Y-2);
    if strcmpi(p.axis,'ij')
        quiver(X, Y, max(Gy(X(1,:),Y(:,1),:),[],3), max(Gx(X(1,:),Y(:,1),:),[],3));
    else
        quiver(X, Y, Gx(X(1,:),Y(:,1)), -Gy(X(1,:),Y(:,1)));
    end
    hold off; title('gradient field')
end

end % end of grdmask
