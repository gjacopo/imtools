%% GRD2GST - Gradient Structure Tensor from directional derivatives.
%
%% Description
% Compute the entries of a (symmetric) gradient structure tensor (GST, aka
% second moment matrix, scatter matrix, Forstner interest operator) of an
% image from its directional derivatives:
%
%       GST = conv ( G(rho), grad(I) * grad(I)^T )
%
% where grad(I) is given by [gx,gy] and G(rho) is a (possibly anisotropic)
% Gaussian kernel of scale rho.
% 
%% Syntax
%       [gx2, gy2, gxy] = GRD2GST(gx, gy);
%       [gx2, gy2, gxy] = GRD2GST(gx, gy, rho, int);
%       [gx2, gy2, gxy, mag, or] = ...
%                 GRD2GST(gx, gy, rho, int, 'Property',propertyvalue, ...);
% 
%% Inputs
% *|I|* : an image with size |(X,Y,C)|, where |C>1| when |I| is multichannel.
%
% *|gx, gy|* : directional derivatives; depending on parameter |axis| (see
%      below), |gx| and |gy| are either:
%
% * the derivatives in I-(vertical oriented NS) and J-(horizontal oriented 
%          OE), like matlab coordinates, when |axis='ij'| (default), or
% * the derivatives in X-(horizontal oriented OE) and Y-(vertical oriented 
%          SN) directions when |axis='xy'|; 
% 
% typically, derivatives can be estimated using |[gy,gx]=GRADIENT(I)| or 
%      |[gx,gy]=GRDSMOOTH(I,...,'axis','ij')|.
%
% *|rho|* : post-smoothing standard deviation; this parameter sets the
%     integration scale for spatial averaging, that controls the size of 
%     the neigxbourhood in which an orientation is dominant; it is used for 
%     averaging the partial directional derivatives of the tensor with a
%     Gaussian kernel; if |rho<0.05|, then no smoothing is performed; default: 
%     |rho=1|.
%
% *|int|* : optional string setting the method used for integrating (spatially
%     averaging) the GST; smoothing is applied, which can be either isotropic,
%     performed in local isotropic neighbourhoods, by setting it to:
%
% * |'fast'| for the (fast) 2D Gaussian smoothing with |IMGAUSSIAN|,
% * |'conv'| for the 2D Gaussian filtering with |GAUSSKERNEL| and |CONVOLUTION|, 
% * |'matlab'| for a Matlab-like implementation of the 2D smoothing with
%          |FSPECIAL| and |IMFILTER|,
%
% or anisotropic, to better capture edges anisotropy, by setting it to:
% 
% * |'ani'| for anisotropic Gaussian filtering along the edges with
%          |HOURGLASSKERNEL| for defining hourglass shaped Gaussian kernels;
%
% if no smoothing needs to be performed, |int| can be set to |false|: this
%     is equivalent to setting |rho=0|; default: |int='fast'| or, equivalently,
%     |int=true|; see function |SMOOTHFILT|.
% 
%% Property [propertyname  propertyvalues]
% *|'axis'|* : string indicating the direction/orientation of the input
%     directional derivatives (see above); in particular, calls to
%     |GRD2GST(gx,gy,...,'axis','ij')| and |GRD2GST(gy,-gx,...,'axis','xy')| 
%     are equivalent; default: |axis='ij'|, ie. the vector orthogonal to the
%     (real) gradient is passed.
%
% *|'hsize'|* : optional filter size; it is estimated depending on |rho|, 
%     typically |hsize=6*rho+1|; default: |hsize=[]|, calculated later on.
%   
% *|'samp'|* : if the gradients are interpolated to avoid aliasing the filters,
%     need to be adapted using this sampling factor; default: |samp=1|.
%
% *|'eign'|* : in the case the tensor norm estimated from the eigenvalues (l1
%     and l2, with l1>l2) is returned (|nargin>=3|), the string |eign| defines
%     the method used for its approximation; it is either (see |GSTFEATURE|):  
%     |'abs'|, |'l1'| (or |'zen'|), |'sum'| (or |'sap'|), |'dif'| (or |'koe'|)
%     or |'ndi'|; default: |eign='l1'|. 
%
% *|'thez', 'sigt'|* : parameters used by the anisotropic filtering with
%     hourglass filters (see |HOUGLASSKERNEL|); default: |thez=8|, |sigt=.4|.
%
% *|'tn'|* : optional flag (|true| or |false|) to normalize the tensor prior
%     to its smoothing; default: |tn=false|.
% 
%% Outputs
% *|gx2, gy2, gxy|* : matrices storing the entries of the GST: 
%
%                   | gx2   gxy |
%           GST  =  |           |     
%                   | gxy   gy2 |
%
% *|mag|* : optional output providing the norm of the tensor.
%
% *|or|* : optional output storing the orientation of the tensor; the double
%      angle representation is used, therefore or in $[-\pi/2,\pi/2]$.
%
%% References
% [Zenzo86]  S. Di Zenzo: "A note on the gradient of a multi-image",
%      CVGIP, 33:116-125, 1986.
%      <http://www.sciencedirect.com/science/article/pii/0734189X86902239>
%
% [Cuma91]  A. Cumani: "Edge detection in multispectral images", CVGIP:
%      Graphical Models and Image Processing, 53(1):40-51, 1991.  
%      <http://www.sciencedirect.com/science/article/pii/104996529190018F>
%
% [VV95]   L. van Vliet and P. Verbeek: "Estimators for orientation and 
%      anisotropy in digitized images", Proc. ASCI, pp. 442-450, 1995. 
%      <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.104.8552>
%
% [Kosch95]  A. Koschan: "A comparative study on color edge detection",
%      Proc. ACCV, pp. 574-578, 1995.
%      <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.31.2648>
%
% [Cuma98]  A. Cumani: "Efficient contour extraction in color images", 
%      Proc. ACCV, pp. 582-589, LNCS 1351, 1998.
%      <http://www.springerlink.com/content/c1750730487puqm0/>
%
% [TD01]  D. Tschumperle and R. Deriche: "Constrained and unconstrained
%      PDE's for vector image restoration", Proc. SCIA, pp. 153-160, 2001.
%      <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.22.4098>
%
% [TD02]  D. Tschumperle and R. Deriche: "Diffusion PDE's on vector-valued
%      images", IEEE Signal Processing Magazine, 19(5):16-25, 2002. 
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1028349>
%
% [Scheun03]  P. Scheunders: "A wavelet representation of multispectral
%      images", Frontiers of Remote Sensing Information Processing, pp.
%      197-224. World Scientific, 2003.
%      <http://ebooks.worldscinet.com/ISBN/9789812796752/9789812796752_0009.html>
%
% [Tschum06]  D. Tschumperle: "Fast anisotropic smoothing of multivalued
%      images using curvature-preserving PDE's", International Journal of
%      Computer Vision, 68(1):65-82, 2006.
%      <http://www.springerlink.com/content/5217329131181uh4/>
%
%% See also  
% Related:
% <GRDSMOOTH.html |GRDSMOOTH|>,
% <GSTSMOOTH.html |GSTSMOOTH|>,
% <../../filter/html/SMOOTHFILT.html |SMOOTHFILT|>.
% Called:
% <GRD2GST_BASE.html |GRD2GST_BASE|>.

%% Function implementation
function [gx2, gy2, gxy, varargout] = grd2gst(gx, gy, varargin)

%%
% parsing and checking parameters

error(nargchk(2, 26, nargin, 'struct'));
error(nargoutchk(3, 5, nargout, 'struct'));

if ~(isnumeric(gx) && isnumeric(gy))
    error('grd2gst:inputerror','matrices are required in input');
end

p = createParser('GRD2GST');   % create an instance of the inputParser class.
% optional parameters
p.addOptional('rho',1., @(x)x>=0); % just for testing
p.addOptional('int', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','conv','fast','ani'}))));
p.addParamValue('hsize',[], @(x)isscalar(x) || isempty(x));
p.addParamValue('samp',1, @(x)isscalar(x) && round(x)==x && x>=1);
p.addParamValue('tn', false, @(x)islogical(x));
p.addParamValue('thez', 8, @(x)isscalar(x) && round(x)==x);
p.addParamValue('sigt', .4, @(x)isscalar(x) && isfloat(x) && x>0);
p.addParamValue('axis','ij',@(x)ischar(x) && any(strcmpi(x,{'ij','xy'})));
p.addParamValue('eign','l1',@(x)ischar(x) && ...
    any(strcmpi(x,{'abs','zen','l1','sap','sum','ndi','dif','koe'})));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%% 
% checking parameters and setting variable

if size(gx) ~= size(gy)
    error('grd2gst:inputerror','matrices must have same dimensions');
end

if strcmp(p.axis,'ij')
    % take the vector orthogonal to the input gradient
    tmp = gx; gx = gy; gy = -tmp;
    % otherwise: horizontal OE and vertical SN derivatives have been passed
    % in gx and gy resp.
end

if islogical(p.int) && p.int,  p.int = 'fast';  end % reset to default

%% 
% main computation

if nargout == 3
    [gx2, gy2, gxy] = grd2gst_base(gx, gy, p.rho, p.int, p.hsize, ...
        p.samp, p.tn, p.thez, p.sigt, p.eign);

else    
    [gx2, gy2, gxy, mag, or] = grd2gst_base(gx, gy, p.rho, p.int, p.hsize, ...
        p.samp, p.tn, p.thez, p.sigt, p.eign);
    varargout{1} = mag;
    if nargout == 5, varargout{2} = or;  end;    
end

%%
% display
if p.disp
    figure,
    ncols = 3; if nargout==3,  nrows = 1; else nrows = 2;  end
    subplot(nrows,ncols,1), imagesc(rescale(gx2,0,1)), axis image off,
    title(['g_x^2 - axis ''' p.axis '''']); colormap gray;
    subplot(nrows,ncols,2), imagesc(rescale(gy2,0,1)), axis image off,
    title(['g_y^2 - axis ''' p.axis '''']); colormap gray;
    subplot(nrows,ncols,3), imagesc(rescale(gxy,0,1)), axis image off,
    title(['g_x \cdot g_y - axis ''' p.axis '''']); colormap gray;
    if nargout>=4
        subplot(nrows,ncols,4), imagesc(rescale(mag,0,1)), axis image off,
        title('mag'), colormap gray
        subplot(2,3,5), imagesc(rescale(or,0,1)), axis image off,
        title('or'), colormap gray
    end
end

end % end of grd2gst

