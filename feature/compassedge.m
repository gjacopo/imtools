%% COMPASSEDGE - Compass edge detection.
%
%% Description
% Nothing else than a wrapper for the Generalized Compass Operator for images
% [RT99a,RT99b,RT01], which contains various information related to edges
% and corners in an image.
% See source code [COMPASS] by Ruzon & Tomasi.
%
%% Algorithm
% The compass operator uses a circular window centered at a junction where
% 4 pixel squares meet. The needle is a diameter at a given orientation. The
% color distributions of the two semicircles are computed, and the distance
% between them is measured using the Earth Mover's Distance (EMD). The 
% maximum EMD over all orientations gives the edge strength and orientation
% at that point. 
% 
%% Syntax
%     [S, O, A, U] = COMPASSEDGE(I);
%     [S, O, A, U] = COMPASSEDGE(I, sigma, 'Property', propertyvalue, ...);
%
%% Inputs
%      
% *|I|* : input image, possibly multichannel, with size |(X,Y,C)|; when the 
%     variable |'gray'| (see below) is set to true, then the function
%     |GREYCOMPASS_MEX| is applied over the different channels of the input
%     image (even in the case |C=3| channels); when |'gray'| is set to false
%     and |C=3| channels, the function |COMPASS_MEX| is used: if |I| is of class
%     |uint8|, the values are treated as RGB values, if |I| is of class double,
%     the values are assumed to be in the CIE-Lab color space. 
%      
% *|sigma|* : (optional) scalar defining the standard deviation of the Gaussian
%     used to weight pixels in the Compass Operator (operator's radius is
%     |3*R|); default: |sigma=1|. 
%
%% Property [propertyname  propertyvalues]
% *|'gray'|* : (optional) logical value set to true when the channels of the
%     image |I| are treated independently using |GREYCOMPASS_MEX|, even in the
%     case |C=3|; default: |gray=false|, ie. if |C=3|, the function |COMPASS_MEX|
%     is applied, |GREYCOMPASS_MEX| is still applied for all other |C| values.
%      
% *|'angles'|* : (optional) angles subtended by the Generalized Compass
%     Operator; if |angles| is a vector, the Compass Operator will detect
%     corners and edges at different angles, and the output arguments will
%     have multiple rows; default: |angles=180| (edges only).
%      
% *|'nwedges'|* : (optional) number of wedges (orientations) in one-quarter
%     of the circle; |nwedges| defines the set of possible orientations and
%     angles, so the values in the |angles| vector must match; default: 
%     |nwedges=6| (|90/6 = 15| degree increments).
%
%% Outputs
% *|S|* : strength of the Compass Operator, lying in the interval |[0,1]|.
%      
% *|O|* : orientation of the edge or corner, in degrees;  for edge detection,
%     this number lies in |[0,180)| and is measured with respect to the
%     positive X-axis; for corners, the interval is |[0,360)| and refers to
%     the orientation of the "right" side of the corner.
%      
% *|A|* : the Abnormality of the Compass Operator at that point; this is the
%     minimum value over all orientations.
%      
% *|U|* : the Uncertainty of the Compass Operator; this is the number of
%     degrees over which, at each orientation, the compass operator gives
%     similar responses.
%
%% References
% [RT99a]  M. Ruzon and C. Tomasi: "Color edge detection with the Compass 
%      Operator", Proc. IEEE CVPR, vol. 2, pp. 160-166, 1999. 
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=784624>
%      
% [RT99b]  M. Ruzon and C. Tomasi: "Corner detection in textured color 
%      images", Proc. ICCV, vol. 2, pp. 1039-1045, 1999.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=790384>
%
% [RT01]  M. Ruzon and C. Tomasi: "Edge, junction, and corner detection
%      using color distributions",  IEEE Transactions on Pattern Analysis
%      and Machine Intelligence, 23(11):1281-1295, 2001. 
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=969118>
% 
% [COMPASS] code available at <http://ai.stanford.edu/~ruzon/compass/>.
%
%% Acknowledgment
% This function calls the mex functions implemented in the original [COMPASS]
% code by Ruzon & Tomasi.
%
%% See also
% Related:
% <EDGECORNER.html |EDGECORNER|>,
% <CANNYEDGE.html |CANNYEDGE|>,
% <CANNYEDGEPROD.html |CANNYEDGEPROD|>,
% <ROTHWELLEDGE.html |ROTHWELLEDGE|>,
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>,
% <SDGDEDGE.html |SDGDEDGE|>,
% <ANISOEDGE.html |ANISOEDGE|>,
% <ELDERZUCKEREDGE.html |ELDERZUCKEREDGE|>,
% <KOETHEDGE.html |KOETHEDGE|>,
% <PETROUEDGE.html |PETROUEDGE|>.
% Called:
% <COMPASSEDGE_BASE.html |COMPASSEDGE_BASE|>.

%% Function implementation
function [S, O, A, U] = compassedge(I, varargin)

%%
% parsing parameters

error(nargchk(1, 16, nargin, 'struct'));
error(nargoutchk(1, 4, nargout, 'struct'));

if ~isnumeric(I)
    error('cannyedge:inputparameter','a matrix is required in input'); 
end

p = createParser('COMPASSEDGE');   
p.addOptional('sigma', 1, @(x)isscalar(x) && x>=0.05);
p.addParamValue('angles', 180, @(x) (isscalar(x) || isvector(x)) && ...
    all(x>=0 & x<=360));
p.addParamValue('nwedges', 6, @(x)isscalar(x) & x>=1);
p.addParamValue('gray', false, @(x)islogical(x));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% calculation

[S, O, A, U] = compassedge_base(I, p.sigma, p.angles, p.nwedges, p.gray);

%%
% display

if p.disp
    figure, imagesc(S), axis image off, title('compass strenght');
    if size(S,3)==1;  colormap gray;  end
end

end % end of compassedge
