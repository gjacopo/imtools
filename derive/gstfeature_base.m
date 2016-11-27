%% GSTFEATURE_BASE - Base function for GSTFEATURE.
% 
%% Syntax
%      [f1, f2, ...] = GSTFEATURE_BASE(gx2, gy2, gxy, lfeature, eign, ex, ey);
%
%% See also
% Related:
% <GSTFEATURE.html |GSTFEATURE|>.
% Called:
% <GSTDECOMP.html |GSTDECOMP|>.

%% Function implementation
%--------------------------------------------------------------------------
function varargout = gstfeature_base(gx2, gy2, gxy, lfeat, eign, ex, ey)

%% 
% internal variables and further testing

% create the list of feature names
if ischar(lfeat),  lfeat = {lfeat};  end
nfeat = numel(lfeat);

%% 
% main computation

if any(strcmp('eigenorm',lfeat)) || any(strcmp('norm',lfeat)) || ...
        any(strcmp('coherence',lfeat)) || any(strcmp('coher',lfeat)) || ...
        any(strcmp('vectorial',lfeat)) || any(strcmp('orvec',lfeat))
    [l1,l2,~,e2] = gstdecomp(gx2, gy2, gxy);         
end


for i=1:nfeat
    switch strtrim(lfeat{i})
        case {'norm','eigenorm'}
            varargout{i} = gstnorm(l1,l2,eign);
        case {'frob','frobenius'}
            varargout{i} = gstfrob(gx2, gy2, gxy);
        case {'orient','orientation'}
            varargout{i} = gstorient(gx2, gy2, gxy);
        case {'orvec','vectorial'}
            theta = gstorient(gx2, gy2, gxy);
            varargout{i} = gstreorient(theta,e2,ex,ey);
        case {'ordir','direction'}
            varargout{i} = gstorientdir(gx2, gy2, gxy);
        case {'inert','inertia'}
            varargout{i} = gstinertia(gx2, gy2, gxy);
        otherwise % case {'coher','coherence'}
            varargout{i} = gstcoherence(l1,l2);
    end
end

end % end of gstfeature


%% Subfunctions

%%
% |GSTNORM| - compute the representative value for the norm.
%--------------------------------------------------------------------------
function tn = gstnorm(l1, l2, eigenmethod)

% initiazilation
if any(strcmp(eigenmethod,{'ni','ndi'}))
    tn = zeros(l1);
    ir = (l1+l2~=0);
    l1 = l1(ir); l2 = l2(ir);
end

switch eigenmethod
    case 'abs'
        tn = abs(l1);
    case {'zen','l1'}
        tn = sqrt(l1);    
    case {'sap','sum'}
        tn = sqrt(l1+l2); 
    case {'dif','koe'} % difference
        tn = sqrt(l1-l2);
    case 'ndi' % normalized difference
        tn(ir) = (l1 - l2) ./ sqrt(l1+l2);
    case 'ni'
        tn(ir) = l1 ./ sqrt(l1+l2);
end

end

        
%%
% |GSTORIENT| - Compute the orientation of the tensor using double angle 
% representation of Wilson's approach.
% The direction of the first eigenvector lambda1 indicates the prominent 
% local orientation, which is equal to the orientation in the image with
% maximum change.
%
%    A. Cumani: "Edge detection in multispectral images" - EQ.(5)  
%    S. di Zenzo: "A note on the gradient of a multi-image"
%    A. Koschan: "A comparative study on color edge detection"
%    L. van Vliet and P. Verbeek: "Estimators for orientation and anisotropy 
%        in digitized images" - EQ.(5)
%
% The orientation of the tensor is given in mathematical positive (counter-
% clockwise) orientation, starting at the (horizontal) X-axis.
%--------------------------------------------------------------------------
function theta = gstorient(gx2, gy2, gxy)

%ir = abs(gx2-gy2) < eps;
%theta = zeros(size(gx2));
%theta(~ir) = .5 * atan2(2*gxy(~ir), gx2(~ir) - gy2(~ir));

theta = .5 * atan2(2*gxy, gx2 - gy2);

% indetermination exists when gxy==0, ie. in :
%   gy==0 & gx~=0  as gx2-gy2 is always >0 
%   gx==0 & gy~=0  as gx2-gy2 is always <0 

% theta is in the range [-pi/2,pi/2]
%theta(ir & gxy>0) = pi/4; 
%theta(ir & gxy<0) = -pi/4;
%theta(ir & gxy==0) = 0;
% all assigments/tests through ir are already performed with atan2

% w is the average value of the magnitude square of the gradient estimate
% w = mean(gx2+gy2);
end


%%
% |GSTORIENTDIR| - Compute the orientation of the tensor.
% Tensor orientation estimation from the main eigenvector coordinates
% recomputed from the tensor partial derivatives 
%
%    A. Cumani: "Efficient contour extraction in color images" - EQ.(5) 
%--------------------------------------------------------------------------
function theta = gstorientdir(gx2, gy2, gxy)

v1 = gx2 - gy2;
v2 = 2 * gxy;
v12 = v1*v1 + v2*v2;

v1 = v1 ./ (sqrt(v12)+eps);  % avoid dividing by 0
theta = atan2(sqrtf(0.5 * (1-v1)) * sign(v2), sqrt(0.5 * (1-v1)));
end


%%
% Analytic solution of principal direction
%
%      den = sqrt(gxy.^2 + (gxx - gyy).^2);
%
% Sine and cosine of doubled angles
%
%      sin2theta = gxy/den;           
%      cos2theta = (gxx-gyy)/den;
%
% The minimum inertia can be estimated as the area moment about the 
% orientation axis:
%
%       imin = 0.5( (gy2+gx2) - (gx2-gy2)*cos2theta - gxy*sin2theta);
%
% The maximum inertia can be found on the perpendicular axis:
%
%       imax = gy2 + gx2 - imin;
%
% The reliability measure of orientation data is given by 
%
%       1.0 - min_inertia/max_inertia.  
%
% The reasoning being that if the ratio of the minimum to maximum inertia 
% is close to 1 we have little orientation information.  

%%
% |GSTREORIENT| - Repeals the undermination in the orientation of the tensor
% using Scheunders approach. Consistent orientation is given by the vector
% field [ex,ey].
%--------------------------------------------------------------------------
function theta = gstreorient(theta, e2, ex, ey) 

% change the direction of the gradient vector:
%   - if the angle (orientation) is <pi (lower part of the quadrant), ADD
%     PI,
%   - if the angle (orientation) is >pi (upper part of the quadrant),
%     SUBSTRACT PI
% this depends on the sign of evy

% ir = ex.*cos(theta) + ey.*evy<0;
% theta(ir) = theta(ir) - pi*sign(evy(ir)); 
s = sign(e2(:,:,2)) - (e2(:,:,2)==0);
theta = theta - pi * s .* (ex.*e2(:,:,1) + ey.*e2(:,:,2)<0);
end


%%
% |GSTFROB| - Compute the Frobenius norm, ie the sum of the squares of the 
% tensor entries.
%--------------------------------------------------------------------------
function frob = gstfrob(gx2, gy2, gxy)                                 
frob = sqrt(gx2 .^ 2 + gy2 .^ 2 + 2 * gxy.^2); 
end


%%
% |GSTINTERTIA| - Compute the maximum contrast in the tensor direction 
%--------------------------------------------------------------------------
function inertia = gstinertia(gx2, gy2, gxy) 
den = (gxy * gxy + (gx2-gy2)*(gx2-gy2)+eps);

inertia =   0.5 * (gx2 + gy2 + (gx2-gy2).^2/den + (gxy*gxy)/den); 
end


%%
% |GSTCOHERENCE| - Compute the coherence (aka anisotropy) index.
% The coherence is a value capable of distinguishing between the isotropic
% and uniform cases; it is obtained as a function of the eigenvalues.
% See H. Wang et al.: "Gradient adaptive image restoration and
% enhancement".
%--------------------------------------------------------------------------
function coherence = gstcoherence(l1,l2) 
coherence = (l1 - l2) ./ (l1+l2+eps);
end
