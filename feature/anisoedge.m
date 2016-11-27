%% ANISOEDGE - Anisotropic diffusion edge detector.
%
%% Description
% Implement the anisotropic diffusion edge detector proposed in [BSMH98].
%
%% Syntax
%     edgemap = ANISOEDGE(I);
%     edgemap = ANISOEDGE(I, sigma);
%     [edgemap,mag] = ANISOEDGE(I, sigma, 'Property', propertyvalue, ...);
% 
%% Inputs
% *|I|* : input image.
%      
% *|sigma|* : optional standard deviation of the Gaussian filter used for
%     smoothing of the image; typical range of values is |[0.5,2]|; default:
%     |sigma=1|.
%
%% Property [propertyname  propertyvalues]
% *|'iter'|* : number of iterations used for smoothing; default: |iter=100|.
%      
% *|'max'|* : optional boolean flag for non-maximal suppression; when set to
%     |true|, the detector is then supplemented with non-maximal suppression 
%     to produce single pixel wide edges; default: |max=true|.
%
%% Outputs
% *|edgemap|* : a binary edge map.
%      
% *|mag|* : optional output image for the gradient magnitude.
%
%% References
% [BSMH98] M. Black, G. Sapiro, D. Marimont and D. Heeger: "Robust 
%      anisotropic diffusion", IEEE Trans. Image Processing, 7(3,):421-432, 
%      1998.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=661192>
%
% [HSSB98]  M. Heath, S. Sarkar, T. Sanocki, and K. Bowyer: "Comparison of
%     edge detectors: a methodology and initial study",  Computer Vision and
%     Image Understanding, 69(1):38-54, 1998.
%     <http://www.sciencedirect.com/science/article/pii/S1077314297905877>
%      
%% Acknowledgment
% This function uses the C functions developped in [HSSB98] by
%     <mailto:kranenbu@bigpine.csee.usf.edu Christine Kranenburg>
% available at <http://marathon.csee.usf.edu/edge/>.
%
%% See also
% Related:
% <EDGECORNER.html |EDGECORNER|>,
% <CANNYEDGE.html |CANNYEDGE|>,
% <CANNYEDGEPROD.html |CANNYEDGEPROD|>,
% <SDGDEDGE.html |SDGDEDGE|>,
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>,
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <ANISOEDGE.html |ANISOEDGE|>,
% <ELDERZUCKEREDGE.html |ELDERZUCKEREDGE|>,
% <ROTHWELLEDGE.html |ROTHWELLEDGE|>,
% <KOETHEDGE.html |KOETHEDGE|>,
% <PETROUEDGE.html |PETROUEDGE|>.
% Called:
% <ANISOEDGE_BASE.html |ANISOEDGE_BASE|>.

%% Function implementation
function [edgemap,varargout] = anisoedge(I, varargin)

%%
% check if possible
if ~exist('anisoedge_mex','file')
    error('anisoedge:mexfile','mex file anisoedge_mex not found');
end

error(nargchk(1, 14, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

if ~isnumeric(I)
    error('anisoedge:inputparameter','matrix required in input'); 
end

%% 
% parsing parameters

p = createParser('ANISOEDGE');   
p.addOptional('sigma',1., @(x)isscalar(x) && x>=0.1); 
p.addParamValue('iter', 100, @(x)isscalar(x) && x>=1);
p.addParamValue('max', true, @(x)islogical(x));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%% 
% main computation

[edgemap,tmpmag] = anisoedge_base(I, p.sigma, p.iter, p.max);

if nargout==2,  varargout{1} = tmpmag;  end;

%%
% display

if p.disp
    figure, imagesc(edgemap), axis image off, title('Anisotropic edge map');
    if size(edgemap,3) == 1, colormap gray;  end;
end

end % end of anisoedge
