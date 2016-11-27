%% FASTCPDA - Fast corner detection.
%
%% Description
% Fast corner detector based on the chord-to-point distance accumulation
% technique of [AL08,ALFR09], itself derived from the curvature-based approach
% introduced in [HY04,HY08].
%
%% Algorithm
% The original algorithm by [HY04,HY08] works by the following steps:
%
% # extract the edge contours from the edge-map, fill the gaps in the
% contours.
% # compute curvature at a low scale for each contour to retain all
% true corners.
% # all of the curvature local maxima are considered as corner candidates,
% then rounded corners and false corners due to boundary noise and details
%  were eliminated.
% # end points of line mode curve were added as corner, if they are not
% close to the above detected corners.
%
% Improvements due to CPDA (Chord-to-Point Distance Accumulation) rely on
% the observation that setting different sigma values based on the curve 
% length is rather impractical, since a different set of edges may be 
% extracted from the test images and the extracted edges often vary in 
% lengths. Therefore, sigma is set to sigma=3 for all curves for the CPDA 
% detector. This is different from [AG08] where sigma was set one of three 
% values 1, 2 and 3 based on the curve-length.
%
%% Syntax
%     [ptcorner,cornermap] = FASTCPDA(I);
%     [ptcorner,cornermap] = FASTCPDA(I, 'Property', propertyvalue, ...);
%     [ptcorner,cornermap] = FASTCPDA(emap);
%     [ptcorner,cornermap] = FASTCPDA(emap, 'Property', propertyvalue, ...);
% 
%% Input
% *|V|* : input image (see also variables |'der'| and |'sig'| below) to analyze
%      or (binary/logical) map of already detected edges.
% 
%% Property [propertyname  propertyvalues]
% *|'gap'|* : maximum gap between the successive edge-points, used to fill 
%     the gaps in the contours; if the gap between the sucessive points is 
%     more than Gap_size, then they would be in two different edges; default: 
%     |gap=1| pixel.
%      
% *|'end'|* : flag to indicate whether the end-points of a curve are detected 
%     as corners (|true|) or not; default: |end=false|.
%      
% *|'thang'|* : threshold on angle orientation; it denotes the maximum obtuse
%     angle that a corner can have when it is detected as a true corner;
%     default: |thang=157|;
%      
% *|'der'|* : in the case, the input image is to be processed, this optional
%     string stores the method used for calculating the edge map: |'matlab'|, 
%     |'vista'| or |'kovesi'|; default: |der='vista'|.
%      
% *|'sig'|* : similarly, this (optional) variable gives standard deviation of
%     the Gaussian filter used for smoothing the image and computing the 
%     curvature; default: |sig=1|.
%      
% *|'reduce'|* : logical value or string defining the way the different channels
%       of a multispectral image are combined into the output edge map:
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
% *|cout|* : |(n,2)| matrix containing the coordinates of the detected corners,
%       where |n| is the number of detected corners.
%      
% *|cd|* : cpda curvature values of the detected corners.
% 
%% References
% [MS98]  F. Mokhtarian and R. Suomela: "Robust image corner detection through 
%      curvature scale space", IEEE Trans. on Pattern Analysis and Machine
%      Intelligence, 20(12):1376-1381, 1998
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=735812>
%
% [HY04]  X.C. He and N.H.C. Yung: "Curvature scale space corner detector
%      with adaptive threshold and dynamic region of support", Proc. ICPR,
%      vol. 2, pp. 791-794, 2004.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1334377>
%      
% [HY08]  X.C. He and N.H.C. Yung: "Corner detector based on global and 
%      local curvature properties", Optical Engineering, 47(5):057008, 2008.
%      <http://spiedigitallibrary.org/oe/resource/1/opegar/v47/i5/p057008_s1>
%
% [AL08]  M. Awrangjeb and G. Lu: "Robust image corner detection based 
%      on the chord-to-point distance accumulation technique", IEEE Trans.
%      on Multimedia, 10(6):1059-1072, 2008.
%      <http://www.gscit.monash.edu.au/gscitweb/loid.php?loid=905299&mimetype=application/pdf>
%      
% [ALFR09]  M. Awrangjeb, G. Lu, C.S. Fraser and M. Ravanbakhsh: "“A Fast
%      Corner Detector Based on the Chord-to-Point Distance Accumulation 
%      Technique", Proc. DICTA, pp. 519-525, 2009.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5384897>
%
%% Acknowledgement
% This is an adaptation of the original code of [AL08,ALFR09] available at:
%  <http://www.mathworks.com/matlabcentral/fileexchange/28207-a-fast-corner-detector-based-on-the-chord-to-point-distance-accumulation-technique>
%
% Part of this code was already from the source code of [HY04,HY08] available at:
%  <http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=7652&objectType=file>
%
%% See also
% Related:
% <CORNER.html |CORNER|>,
% <EDGECORNER.html |EDGECORNER|>,
% <SUSANCORNER.html |SUSANCORNER|>,
% <HARRISCORNER.html |HARRISCORNER|>,
% <FASTCORNER.html |FASTCORNER|>,
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>.
% Called: 
% <FASTCPDA_BASE.html |FASTCPDA_BASE|>,
% <CANNYEDGE_BASE.html |CANNYEDGE_BASE|>.

%% Function implementation
function [ptcorner,cornermap] = fastcpda(V, varargin)

%%
% parsing and checking parameters

error(nargchk(1, 21, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

% mandatory parameter
if ~(isnumeric(V) || islogical(V))
    error('fastcpda:inputerror','numeric or logical matrix required in input'); 
end

p = createParser('FASTCPDA');   % create an instance of the inputParser class.
% optional parameters
p.addParamValue('thang',157, @(x)x>=0 && x<=180); 
p.addParamValue('gap',1, @(x)isscalar(x) && x>=1);
p.addParamValue('end',false,@(x)islogical(x));
% used in case an edge map has to be produced first
p.addParamValue('sig', 1, @(x)isscalar(x) && x>=0.05);
p.addParamValue('der', 'vista', @(x)ischar(x) && ...
    any(strcmpi(x,{'matlab','vista','kovesi'})));
p.addParamValue('reduce', false, @(x)islogical(x) || ...
    (ischar(x) && any(strcmpi(x,{'igray','imax','isum','gmax','eor'}))));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% calculation

if ~islogical(V)
    V = cannyedge_base(V, p.sig, p.der, [], p.reduce);
end

[ptcorner, cornermap] = fastcpda_base(V, p.gap, p.thang, p.end);

%%
% display

if p.disp
    tmp = V; % to show corners on the input image
    for c=1:size(ptcorner,1)
        for i=1:size(ptcorner{c},1)
            tmp(:,:,c) = mark(tmp(:,:,c),ptcorner{c}(i,1),ptcorner{c}(i,2),3);
        end
    end    
    figure, imagesc(tmp), axis image off, title('fast CPDA corner map');
    if size(tmp,3) == 1,   colormap gray;  end
end


end % end of fastcpda


%% Subfunction

%%
% |MARK| - Show corners into the output images or into the edge-image.
%--------------------------------------------------------------------------
function I1 = mark(I, x, y, w)
[M,N,~] = size(I);                                                  

I1 = I;

if isa(I,'logical')
    I1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:) = ...
        (I1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)<1);
    I1(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:) = ...
        I(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:);

else
    I1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:) = ...
        (I1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)<128)*255;
    I1(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:) = ...
        I(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:);
end

end % end of mark
