%% CORNER - Corner detection techniques.
%
%% Description
% Performs the implementation of and the call to various corner detection methods.
%
%% Syntax
%     [cornermap, ptcorner] = CORNER(I, method, 'Property', propertyvalue, ...);
%
%% Inputs
% *|I|* : input image, possibly multichannel.
%      
% *|method|* : (optional) string defining the method used for corner extraction;
%      it is either:
%      
% * |'harris'| or |'noble'| for Harris-Forster and Noble detections resp.
%          using the function |HARRISCORNER|,
% * |'susan'| for Susan detection [SB97] with function |SUSANCORNER|,
% * |'fast9'|, |'fast10'|, |'fast11'|, |'fast12'|, for any FAST corner detection
%          using the function |FASTCORNER|,
% * |'cpda'| for curvature-based detection [ALFR09] using the function
%          |FASTCPDA|;
%      
% note: currently, curve corner detection (|'cur'|, outperformed by |'cpda'|), 
% compass corner detection (|'com'|, very slow) are not included in this function,
% but could easily be added.
%
%% Property [propertyname  propertyvalues]
% *|'thres'|* : threshold parameter used by most of the corner detectors; see
%      respective calls; default: |thres=0.001| (as used by Harris detector).
%      
% *|'reduce'|* : logical value defining the way the different channels
%      of a multispectral image are combined into the corner map; when set 
%      to |true|, the detector is applied channel by channel and the output
%      corner map is taken as the logical |OR| of the output corner maps of
%      the different channels; default: |reduce=true|.
%      
% *|'gap', 'thang'|* : parameters used by |'cpda'|, maximum gap between the
%      successive edge-points and threshold on angle orientation;; default:
%      |gap=1| pixel and |thang=157| degrees.
%      
% *|'kappa', 'radius'|* : parameters used by |'har'| and/or |'nob'|; default: 
%      |kappa=0.06| and |radius=3|.
%      
% *|'sig'|* : numeric variable used when the selected method is any of |'cpda'|,
%      |'har'| or |'nob'| as the standard deviation for smoothing the 1st order
%      directional derivatives prior to edge detection; default: |sig=1|.
%      
% *|'rho'|* : numeric variable used when the selected method is any of |'har'|,
%      or |'nob'| as the standard deviation for averarging the 2nd order
%      directional derivatives prior to edge detection; default: |rho=1|.
%
%% Outputs
% *|cornermap|* : logical 2D map of corners.
%      
% *|ptcorner|* : matrix storing corners' coordinates.
%
%% References
%   [HS88]  C.G. Harris and M.J. Stephens: "A combined corner and edge 
%      detector", Proc. Vision Conference, pp 147-151, 1988.
%      <http://www.bmva.org/bmvc/1988/avc-88-023.pdf>
%      
%   [Noble89]  A. Noble: "Descriptions of Image Surfaces", PhD thesis, p45, 
%      Oxford University 1989.
%      
%   [SB97]  S.M. Smith and J.M. Brady: "{SUSAN} - A new approach to low  
%      level image processing", International Journal of Computer Vision,
%      23(1):45-78, 1997.
%      <http://perception.inrialpes.fr/Publications/2000/SMB00/Schmid-ijcv00.pdf>
%      
%   [HY04]  X.C. He and N.H.C. Yung: "Curvature scale space corner detector
%      with adaptive threshold and dynamic region of support", Proc. ICPR,
%      vol. 2, pp. 791-794, 2004.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1334377>
%      
%   [HY08]	X.C. He and N.H.C. Yung: "Corner detector based on global and 
%      local curvature properties", Optical Engineering, 47(5):057008, 2008.
%      <http://spiedigitallibrary.org/oe/resource/1/opegar/v47/i5/p057008_s1>
%      
%   [AL08]  M. Awrangjeb and G. Lu: "Robust image corner detection based 
%      on the chord-to-point distance accumulation technique", IEEE Trans.
%      on Multimedia, 10(6):1059-1072, 2008.
%      <http://www.gscit.monash.edu.au/gscitweb/loid.php?loid=905299&mimetype=application/pdf>
%      
%   [ALFR09]  M. Awrangjeb, G. Lu, C.S. Fraser and M. Ravanbakhsh: "“A Fast
%      Corner Detector Based on the Chord-to-Point Distance Accumulation 
%      Technique", Proc. DICTA, pp. 519-525, 2009.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5384897>
%
%% See also
% Related:
% <EDGECORNER.html |EDGECORNER|>,
% <SUSANCORNER.html |SUSANCORNER|>,
% <HARRISCORNER.html |HARRISCORNER|>,
% <FASTCORNER.html |FASTCORNER|>,
% <FASTCPDA.html |FASTCPDA|>,
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>.
% Called: 
% <CORNER_BASE.html |CORNER_BASE|>.

%% Function implementation
function [cornermap, ptcorner] = corner(I, varargin)

%% 
% parsing parameters

error(nargchk(1, 28, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(I)
    error('corner:inputerror','matrix required in input'); 
end

% optional parameters
p = createParser('CORNER');  
p.addRequired('method',@(x)ischar(x) || ...
    any(strcmpi(x,{'harris','noble','susan','cpda', ...
    'fast9','fast10','fast11','fast12'})));
% used by FASTCORNER
p.addParamValue('nonmax', true, @(x)islogical(x));
% % used by SUSANCORNER
% p.addParamValue('md', 'flat', @(x)(ischar(x) && strcmp(x,flat)) || ...
%     (isscalar(x) && x>0));
% p.addParamValue('n', false, @(x)isscalar(x) && (x==true||x==false));
% p.addParamValue('q', true, @(x)isscalar(x) && (x==true||x==false));
% used by HARRISCORNER
p.addParamValue('kappa', 0.06, @(x)isscalar(x) && isfloat(x) && x>0);
p.addParamValue('radius', 3, @(x)isscalar(x) && round(x)==x && x>0);
% used by FASTCPDA
p.addParamValue('thang',157, @(x)x>=0 && x<=180); 
p.addParamValue('gap',1, @(x)isscalar(x) && x>=1);
% used by FASTCPDA and HARRISCORNER
p.addParamValue('sig', 1, @(x) isnumeric(x) || ...
    (isscalar(x) && isfloat(x) && x>0));
p.addParamValue('rho', 1, @(x) isnumeric(x) || ...
    (isscalar(x) && isfloat(x) && x>0));
% used by COMPASS
% used by all
p.addParamValue('thres', 0.001, @(x)isscalar(x) && x>=0);
p.addParamValue('reduce',true, @(x)islogical(x));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            


%%
% main calculation

[cornermap, ptcorner] = ...
    corner_base(I, p.method, p.thres, p.kappa, p.radius, ...
    p.nonmax, p.gap, p.thang, p.sig, p.rho, p.reduce );

if p.disp
end


end % end of corner

