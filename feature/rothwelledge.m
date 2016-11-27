%% ROTHWELLEDGE - Subpixel edge detector. 
%
%% Description
% Implementation of the edge detector proposed by Rothwell et al. [RMHN95].
%
%% Syntax
%     edgemap = rothwelledge(I);
%     edgemap = rothwelledge(I,sigma);
%     edgemap = rothwelledge(gx,gy);
%     [edgemap,mag,orient] = ...
%                 rothwelledge(I, sigma, 'low', low, 'a', alpha);
%     [edgemap,mag,orient] =  ...
%                 rothwelledge(gx, gy, 'low', low, 'a', alpha);
% 
%% Inputs
% *|I|* : input image or partial derivative in X-direction, depending on the
%     subsequent option passed.
%      
% *|V|* : optional parameter storing:
%
% * either the standard deviation of the Gaussian filter used for smoothing
%          of the image; typical range of values is |[0.5,2]|;
% * or a matrix with the partial derivative in Y-directions; in that case,
%          edge detection is performed through subpixel interpolation over  
%          the already computed partial derivatives passed in input as
%          |(I,V)|.
% 
% default: |V=1| as the standard deviation (an image is assumed to be passed).
%
%% Property [propertyname  propertyvalues]
% *|'low'|* : unique threshold value used by the Rothwell technique; typical
%     range of values is |[3,18]|; default: |low=5|.
%      
% *|'a'|* : typical range of values is |[0.8,0.95]|; default: |a=0.8|. 
%
% *|'samp'|* : to perform |(* samp)| interpolation of the input image to
%     avoid aliasing; default: |samp=1|.
%      
% *|'reduce'|* : logical value or string defining the way the different channels
%     of a multispectral image are combined into the output edge map (see
%     |ROTHWELLEDGE_BASE|); it can be either:
%      
% * |'igray'| : the input RGB image is converted to a gray image using the
%          function |RGB2GRAY|,
% * |'imax'|, |'isum'| : the input image (any dimension) is converted to a
%          gray image by taking the sum and the max over the different
%          channels resp.,
% * |'gmax'| : gradients are computed for the different channels and their
%          local pixelwise max is given as a single input to the edge
%          detector;
% * |'eor'| : calculations are made like for a multispectral image (as if
%          |reduce=false|), but the final edge map is taken as the logical 
%          |OR| of the output edge maps of the different channels;
%      
% in |true| case, it is set to |'sum'|; default: |reduce=false|, ie. no
%     'combination' is used, a multichannel map is output; this parameter
%     is naturally ignored when |I| is a scalar image.
%
%% Outputs
% *|edgemap|* : a binary edge map.
%      
% *|mag, orient|* : optional output images for the gradient magnitude and the
%     gradient orientation resp.
%
%% Acknowledgment
% This function uses the C function developped Heath et al. for the comparative
% study in [HSSB97]; the mex file in |ROTHWELL_MEX| calls directly the code
% made available by the authors in the page: 
% <ftp://figment.csee.usf.edu/pub/Edge_Comparison/source_code/rothwell.src>
%
%% References
% [RMHN95] C. Rothwell, J. Mundy, B. Hoffman and V.-D. Nguyen: "Driving
%      Vision by Topology", Proc. International Symposium on Computer Vision,
%      pp. 395-400 - Techn. report 2444, INRIA, 1995.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=477034&tag=1>
%
% [HSSB97]  M. Heath, S. Sarkar, T. Sanocki, T. and K.W. Bowyer: "A robust
%      visual method for assessing the relative performance of edge-detection
%      algorithms", IEEE Trans. on Pattern Analysis and Machine Intelligence,
%      19(12):1338-1359, 1997.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=643893>
%
%% See also
% Related:
% <EDGECORNER.html |EDGECORNER|>,
% <CANNYEDGE.html |CANNYEDGE|>,
% <CANNYEDGEPROD.html |CANNYEDGEPROD|>,
% <CANNYEDGEMAP.html |CANNYEDGEMAP|>,
% <SDGDEDGE.html |SDGDEDGE|>,
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>,
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <ANISOEDGE.html |ANISOEDGE|>,
% <ELDERZUCKEREDGE.html |ELDERZUCKEREDGE|>,
% <KOETHEDGE.html |KOETHEDGE|>,
% <PETROUEDGE.html |PETROUEDGE|>.
% Called:
% <ROTHWELLEDGE_BASE.html |ROTHWELLEDGE_BASE|>.

%% Function implementation
function [edgemap,varargout] = rothwelledge(I, varargin)

%%
% check if possible
if ~exist('rothwelledge_mex','file')
    error('mex file rothwell_mex not found');
end

error(nargchk(1, 18, nargin, 'struct'));
error(nargoutchk(1, 3, nargout, 'struct'));

if ~isnumeric(I)
    error('rothwelledge:inputparameter','a matrix is required in input'); 
end

%% 
% parsing parameters

p = createParser('ROTHWELLEDGE');   
p.addOptional('V', 1., @(x)isnumeric(x) || (isscalar(x) && x>=0.05)); 
p.addParamValue('low', 5, @(x)isscalar(x) && x>=0 && x<=20);
p.addParamValue('a', 0.8, @(x)isscalar(x) && x>=0 && x<=1);
p.addParamValue('samp', 1, @(x)isscalar(x) && round(x)==x && x>=1 && x<=5);
p.addParamValue('reduce', false, @(x)islogical(x) || ...
    (ischar(x) && any(strcmpi(x,{'igray','imax','isum','gmax','eor'}))));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% call the mex file
[edgemap,or,mag] = rothwelledge_base(I, p.V, p.low, p.a, p.samp, p.reduce);

if nargout>=2
    varargout{1} = mag;
    if nargout==3
        varargout{2} = or;
    end
end

if p.disp
    figure, imagesc(edgemap), axis image off, title('Rothwell edge map');
    if size(edgemap,3) == 1, colormap gray;  end;
end

end % end of rothwelledge
