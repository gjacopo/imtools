function [P, N, T] = ...
    im2potential_base(I, method, a, rho, sigma, der, int, samp, eign)
% IM2POTENTIAL_BASE - Given an input image, design a metric to be used as a
% potential (cost) function in front propagation over its domain. Typically, 
% we extract norm and directional information from its derivatives to build
% the metric. In particular, (iso- or anisotropic) metrics can be designed
% from the gradient or the gradient structure tensor of the input image.
%
%   P = IM2POTENTIAL_BASE(I, method);
%   P = IM2POTENTIAL_BASE(I, method, a, rho, sigma, der, int, samp, eign);
% 
%   P = IM2POTENTIAL_BASE(I, 'uni', a);
%   P = IM2POTENTIAL_BASE(I, 'pix'/'pixinv', a);
%   P = IM2POTENTIAL_BASE(I, 'grd*'/'iso', a, sigma);
%
%   P = IM2POTENTIAL_BASE(I, 'gst*'/'iso'/'ani', a, rho, sigma, der, int, samp, eign);
%   P = IM2POTENTIAL_BASE(I, 'gst*'/'iso'/'ani', a, rho, sigma, der, int, samp, N);
%   P = IM2POTENTIAL_BASE(T, 'gst*'/'iso'/'ani', a, eign);
%   P = IM2POTENTIAL_BASE(T, 'gst*'/'iso'/'ani', a, N);   
% 
% Inputs:
%   I : input image of size [X,Y,C], possibly multichannel when C>1.
%   T : input tensor matrix, typically the gradient structure tensor (GST);
%     this avoids redoing calculations inside the function.
%   method : string defining the method used for computing the potential
%     (and the metric) derived from the image; it is either based on the image
%     itself by setting it to:
%        - 'pix', then the potential is scalar and equal everywhere to the
%          image (note that the higher the value at a pixel, the higher its
%          potential, the faster the propagation through it),
%        - 'pixinv', then the potential is scalar and equal to the inverse 
%          of the image (inversely: the higher the value at a pixel, the
%          lower its potential, the slower the propagation through it),
%     or based on the gradient of the image, by setting it to:
%        - 'grd', then the potential is vectorial equal to the image gradient,
%        - 'grdorth' or 'ani', then the potential is vectorial equal to the
%          orthogonal vector to the gradient vector,
%        - 'grdn' (or 'isoinv'), then the potential is scalar and equal to
%          the gradient norm,
%        - 'grdninv' (or 'iso'), then the potential is scalar and equal to
%          the inverse of the gradient norm;
%     whenever the image is a scalar (graylevel) image; when the image is
%     multispectral (C>1), the GST is used to derive the potential function,
%     which will depend on the chosen method [GSD10]:
%        - 'iso' (or 'gstninv') for a scalar potential set to the inverse of
%          the GST norm,
%        - 'gst' for a tensorial potential set to the GST itself (same
%          eigendecomposition),
%        - 'gstorth' for a tensorial potential set to the orthogonal tensor
%          to the GST (same eigenvalues, orthogonal eigenvectors),
%        - 'ani' for a tensorial tensor set to the orthogonal of the GST and
%          scaled by N (scaled eigenvalues, orthogonal eigenvector),
%        - 'gstcoh' for a tensorial tensor set to the orthogonal of the GST
%          and scaled by the coherence (see GSTFEATURE),
%        - 'gstiso' for a unitary tensorial potential set to the normalized
%          GST (same eigenvectors, normalized eigenvalues),
%     where N is a scaling function passed externally or estimated as the
%     GST norm (see below).
%   a : exponent used for amplyfying the strenght of the cost function
%     (cases 'pix', 'pixinv', 'grdn', 'iso') or the strenght of the scaling
%     function (cases ''gstninv', 'gstnorm'); default: a=1.
%   rho : integration scale for computing the GST; default: rho=1.
%   sigma : differentiation scale for estimating the directional derivatives;
%     default: sigma=1
%   der : string defining the method of pre-smoothing/differentiation used
%     for estimating the directional derivatives of the input image; it is 
%     either (see GRDSMOOTH): 'matlab', 'vista', 'fast', 'conv', 'fleck', 
%     'tap5', 'tap7', 'sob', 'opt' or 'ana'; default: der='fast'.
%   int : string defining the method used for the post-smoothing of the GST;
%     it is either (see GRD2GST): 'matlab', 'conv' or 'fast' for isotropic 
%     Gaussian smoothing, or 'ani' for anisotropic Gaussian (using hour-
%     glass shaped Gaussian kernels) along the edges; this latter better
%     captures edges anisotropy; default: int='fast'.
%   samp : scalar used as a sampling rate for the gradient when estimating
%     the GST; default: samp=1.
%   eign : in the case the tensor norm estimated from the eigenvalues (l1
%     and l2, with l1>l2) is to be estimated , the string eign defines the
%     method used for its approximation; it is either (see GSTFEATURE): 'l1' 
%     (or 'zen'), 'abs', 'sum' (or 'sap'), 'dif' (or 'koe') or 'ndi'; default: 
%     eign='l1'. 
%   N : additional scaling (or pilot) function; it is a matrix of size [X,Y]
%     used to scalenormalize the eigenvalues of the GST when the potential
%     is derived from it; it counter-balances the influence of the
%     eigenvalues in the potential definition; typically, it is set to the
%     input image, or the gradient norm; default N = [] and the gradient 
%     norm is estimated and used as the scale function.
% 
% Output:
%   P : potential function derived from the image (see also function 
%     POTENTIAL2FRONT_BASE); it can be:
%        - a vector of size [n x m x2] representing a vector field used to
%          define an anisotropic metric when the image is scalar
%        - a matrix of size [n x  m] representing a scalar potential field;
%          it will be used for defining an isotropic metric,
%        - a tensor matrix of size [n x m] representing a tensor field when
%          the image is multispectral, and used for defining an anisotropic
%          Riemannian metric.
%
% Reference:
%   [PC09]  G. Peyre, and L. Cohen: "Geodesic methods for shape and surface
%      processing", in "Advances in Computational Vision and Medical Image
%      Processing: Methods and Applications", vol. 13 of "Computational 
%      Methods in Applied Sciences", pp. 29-56, Springer, 2009.
%   [GSD10]  J. Grazzini, P. Soille and S. Dillard: "Multichannel image 
%      regularisation using anisotropic geodesic filtering", Proc. ICPR,
%      pp. 2664-2667, 2010.
%
% TODO: implement the loop for local window analysis inside the mex for
% efficiency!
%
% See also 
% Related: IM2FRONT_BASE, POTENTIAL2FRONT_BASE, FMMANISOPROPAGATION_BASE, 
%   FMMISOPROPAGATION_BASE, FMM_BASE --
% Called: GRDSMOOTH_BASE, GSTSMOOTH_BASE, GSTFEATURE_BASE, GSTDECOMP.

% note that not all the cases below are listed in the cases offered to the
% user

%% Check/set parameters

error(nargchk(2, 9, nargin, 'struct'));
error(nargoutchk(1, 3, nargout, 'struct'));

% we accept a variable number of inputs
if nargin<9,     eign = 'l1';
    if nargin<8,     samp = 1;
        if nargin<7,     int = 'fast';
            if nargin<6,     der = 'fast';
                if nargin<5,     sigma = 1;
                    if nargin<4,     rho = 1;
                        if nargin<3,     a = [1 2];  end
                    end
                end
            end
        end
    end
end
% still, in the rest of the code we take into account the number of
% variables that have been entered!

if length(a)==1,  a = [a a];
% elseif length(a)==2 &&  a(1)>a(2)
%     a = [a a(1)];
%     a(1) = []; % invert the order
end

d = nb_dims(I);
C = size(I,3);
cons = 1;  % cons = eps;

% dummy settings in some cases
if nargout>=2,  N = [];
    if nargout==3,  T = [];  end
end


%% Define the potential function

if nargin<=3 && strcmpi(method,'uni')  % uniform isotropic potential
    % no information regarding the image is used here
    % compute the normalized unitary tensor: equivalent to Euclidean
    % isotropic unitary eigenvector.
    P = a(1) * ones(size(I)); % Euclidean when a(1) = 1
    % note: this is equivalent to building a tensor metric with the following
    % eigendecomposition:
    % l1 = ones(size(l1)); l2 = ones(size(l2));
    % e1(:,:,1) = 0; e1(:,:,2) = 1; e2 = cat(3, -e1(:,:,2), e1(:,:,1));

elseif C==1 && any(strcmpi(method,{'pix','pixinv'}))
 
    switch method            
        case 'pix' % pixel based value potential
            % the highest the greylevel, the lowest the potential
            P = (eps + I).^a(1);
            
        case 'pixinv'  % pixel based value potential
            % the lowest the greylevel, the highest the potential
            P = 1 ./ (cons + I).^a(1);
    end

    
elseif C==1 && any(strcmpi(method(1:3),{'grd','iso'}))
   
    if nargin<=3 && any(strcmpi(method,{'grdninv','iso','grdn'}))
        N = I;
    else % note that we use rho here!
        [gx, gy, N] = grdsmooth_base(I, rho, 'fast', [], 'ij'); 
    end 
    
    switch  method
        
        case 'grd'
            % vector field
            P(:,:,1) = gx;  P(:,:,2) = gy;

        case {'grdorth','ani'}
            % orthogonal vector field
            P(:,:,1) = gy;  P(:,:,2) = gx;

        case {'grdn','isoinv'} % gradient attracting potential
            P = (eps + N).^a(1);

        case {'grdninv','iso'} % gradient avoiding potential
            % the propagation will be isotropic and based on a scalar field
            % defined as the inverse of the gradient norm (the higher the
            % gradient N, the lower P, the slower the front will propagate)
            P = 1 ./ (cons + N).^a(1);
            
    end
    
    % note that calling im2potential_base(X,'pixinv',a) or  
    % im2potential_base(X,'grdninv',a) is equivalent, but calling
    % im2potential_base(X,'grdninv',a, rho) is not.
    % similarly, it is equiavelent to call im2potential_base(X,'pix',a)
    % or im2potential_base(X,'grdn',a), but not the call to
    % im2potential_base(X,'grdn',a,rho) 
    
elseif any(strcmpi(method(1:3),{'gst','iso','ani'})) % the input can be any dimension
    % compute the structure tensor

    if d<4
        P = gstsmooth_base(I, rho, sigma, der, int, samp, ...
            [], false, false, 8, .4); % default choices
        if nargout==3,  T = P;  end  % may be useful later in some code...
        
    elseif nargin<=4 && d==4  % we suppose we passed the already estimated tensor
        P = I;
        % trick to pass only 4 parameters when the first one is a tensor
        % and only the norm needs to be computed (because then variables
        % rho, sigma, der, int, samp are not needed) by  passing eign in rho
        % or the norm itself in rho
        if ischar(rho) || (~isscalar(rho) && isequal(size(rho),size(I(:,:,1,1))))
            eign = rho;
        end
    end
    
    if ischar(eign)
        % compute the norm of the structure tensor
        N = gstfeature_base(P(:,:,1,1), P(:,:,2,2), P(:,:,1,2), ...
            'norm', eign, [], []);
        % N a can also be used for tensor normalization
        % N = rescale(N,0,1); % N = rescale(N,0,1-eps);
        
    elseif isequal(size(eign),size(I(:,:,1,1))) % whatever I may be (tensor
        % or image), we assume we passed a pilot image in the variable eign
        N = eign;
    end
    
    if any(strcmpi(method,{'iso','gstninv','gstn'}))
        switch method
            case {'iso','gstninv'}
               % we return the inverse of the tensor norm for isotropic propagation
                % (the higher the gradient norm value of a pixel, the lower its
                % potential, the slower the propagation through it)
                 P = 1 ./ (cons + N).^a(1);

            case 'gstn'
                P = (eps + N).^a(1);
        end
        
        return;
    end
    
    % eigendecomposition of the structure tensor
    [l1, l2, e1, e2] = gstdecomp(P);
    
    % % normalize wrt l1's extrema 
    % m = min(l2(:));   mm = max(l2(:));
    % l1 = (l1 - m) / (mm - m);
    % l2 = (l2 - m) / (mm - m);
 
    switch  method
        
        % case 'gstninv'
        %     % tensorial unitary (eigenvectors) and scaled by N
        %     l1 = 1 ./ (cons + N).^a(1);
        %     l2 = eps;
        %     %e1(:,:,1) = 0; e1(:,:,2) = 1;
        %     e2 = cat(3, -e1(:,:,2), e1(:,:,1));
        % % note that in the case the image is scalar (graylevel, C=1)
        % % and N is computed using GSTFEATURE, the results obtained
        % % with methods 'gdrninv' and 'gstninv' are equivalent
        
        case 'gst'
            % do nothing: we will use the GST itself as the tensor metric
            
        case 'gstorth'
            % define the tensor orhogonal to the GST
            tmp = e1;
            e1 = e2;
            e2 = tmp;
            
        case 'ani'
            % the propagation is anisotropic and based on the normalized
            % gradient structure tensor
            l1 = l1 ./ (cons + N).^a(1);
            l2 = 1 ./ (cons + N).^a(2);
 
        case 'gstn1'
            l1 = l1 ./ (cons + N).^a(1);
            l2 = l2 ./ (cons + N).^a(2);
            
            % 'ani' and 'gstn1' perform better
            
        case 'gstn2' % this is exactly Eqs.(3) and (4) of [GSD10]
            l1 = 1 ./ (cons + N).^a(1);
            l2 = 1 ./ (cons + N).^a(2); 
            
        case 'gstn3'
            tmp = l1;
            l1 = l2 ./ (cons + N).^a(1);
            l2 = tmp ./ (cons + N).^a(2);
            
        case 'gstcoh'
            l1 = (l1 + l2) ./ (l1 - l2 + eps);
            l2 = 1;
                        
        case 'gstiso'
            % conserve the GST direction only, set eigenvalues to 1
            l1 = ones(size(l1));
            l2 = ones(size(l2));
            
    end

    % recompose the structure tensor with modified eigenvalues and/or
    % eigenvectors
    P = gstdecomp(l1, l2, e1, e2);  % new tensor metric
    
else
    error('im2potential_base:methoderror', ...
        ['unknown method ' method ' or incompatible combination of input parameters'])
end

end