%% SDGDEDGE - Edge detection based on gradient second derivative.
%
%% Description
% Edge detection filter with 2D kernel based on the 2nd derivative in the
% gradient direction or its combination with the Laplacian. This function
% convolves the input image with set 2D Gaussian kernels and then combines 
% the results in accordance to the formula in [MH86].
%
%% Syntax
%     edgemap = SDGDEDGE(I);
%     edgemap = SDGDEDGE(I, sigma, 'propertyname',propertyvalues ...);
%
%% Inputs
% *|I|* : input image of size |(X,Y,C)|, possibly multichannel with |C>1|.
%
% *|sigma|* : optional standard deviation of the Gaussian filter used for
%     smoothing of the image; default: |sigma=1|.
%
%% Property [propertyname  propertyvalues]
% *|'plus'|* : optional boolean flag for using the 2nd derivative in the
%     gradient direction combined with the Laplacian filter (filter 'PLUS',
%     see [VV94]) when set to |true|; the regular SDGD is used otherwise;
%     default: |plus=false|. 
%
% *|'thres'|* : threshold (in |[0,1]|) used for estimating the binary map from
%     the derivative output filter; default: |thres=0.1|.
%
% *|'hsize'|* : size of the window of the filtering kernel; default: |hsize|
%     is estimated depending on |sigma|.
%
%% Outputs
% *|edgemap|* : the output edge map.
%
%% Acknowledgment
% <mailto:sergei.koptenko{at}resonantmedical.com  Sergei Koptenko>, Resonant
% Medical, Montreal (Qc., Canada), <www.resonantmedical.com>.  
%
%% References
% [MH86]  D. Marr and E.C. Hildreth: "Theory of edge Detection", Proceedings
%      of the Royal Society, 207:187-217, 1980.
%      <http://rspb.royalsocietypublishing.org/content/207/1167/187>
%
% [VV94]  P.W. Verbeek and L.J. Vliet: "On the location error of curved
%      edges in low-pass filtered 2-D and 3-D images", IEEE Trans. on Pattern 
%      Analysis and Machine Intelligence, 16(7):726-733, 1994.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=297954>
%
%% See also
% Related:
% <EDGECORNER.html |EDGECORNER|>,
% <CANNYEDGE.html |CANNYEDGE|>,
% <CANNYEDGEPROD.html |CANNYEDGEPROD|>,
% <ROTHWELLEDGE.html |ROTHWELLEDGE|>,
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>,
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <ANISOEDGE.html |ANISOEDGE|>,
% <ELDERZUCKEREDGE.html |ELDERZUCKEREDGE|>,
% <KOETHEDGE.html |KOETHEDGE|>,
% <PETROUEDGE.html |PETROUEDGE|>.
% Called:
% <SDGDEDGE_BASE.html |SDGDEDGE_BASE|>.

%% Function implementation
function edgemap = sdgdedge(I, varargin)

%% 
% parsing parameters 

error(nargchk(1, 18, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

if ~isnumeric(I)
    error('sdgdedge:inputparameter','a matrix is required in input'); 
end

p = createParser('SDGDEDGE');   
p.addOptional('sigma', 1., @(x)isscalar(x) && x>=0.05);
p.addParamValue('plus', false, @(x)islogical(x));
%p.addParamValue('method', 'sdgd', @(x)ischar(x) && ...
%    any(strcmpi(x,{'sdgd','plus'})));
p.addParamValue('thres', 0.1, @(x)isscalar(x) && x<=1 && x>0);
p.addParamValue('hsize',[], @isscalar);
p.addParamValue('mu', 0, @isscalar); % not used

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% main computation

edgemap = sdgdedge_base(I, p.plus, p.sigma, p.mu, p.thres, p.hsize);

%%
% display

if p.disp
    figure, imagesc(edgemap), axis image off, title('SDGD edge map');
    if size(edgemap,3) == 1, colormap gray;  end;
end

end
