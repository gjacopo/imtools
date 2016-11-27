%% DIRGAUSSKERNEL - Oriented Gaussian kernel.
%
%% Description
% Compute an oriented Gaussian kernel using the approach of [GSW03] and 
% possibly using the code of [ANIGAUSS].
%
%% Syntax
%       H = DIRGAUSSKERNEL(s1, s2, theta, m, method);
%     
%% Inputs
% *|s1|* : vector of size |(1,n)| (or |(n,1)|) of standard deviation(s) along
%      the first (short) axis.
%
% *|s2|* : vector of standard deviation(s) along the (long) second axis, with
%      same dimension as |s1|.
%
% *|theta|* : vector of orientation(s) of the kernel, with same dimension as
%      |s1| (in radians).
%
% *|m|* : maxium size of the kernel (typ. |m=128|).
%
% *|method|* : string setting the implementation used for computing the
%      anisotropic filters; it is:
% 
% * |'geu'| when using the function |ANIGAUSS| implementing the recursive 
%           anisotropic Gaussian proposed in [GSW03,ANIGAUSS],
% * |'2d'| when using a standard 2D implementation of the Gaussian.
%
%% Output
% *|h|* : (array of) matrix(ces) (all with maximum size |(m,m)|) storing the
%      |n| oriented Gaussian filters.
%  
%% References
% [GSW03]  J.M. Geusebroek, A.W.M. Smeulders, and J. Van de Weijer: "Fast
%      anisotropic gauss filtering", IEEE Trans. on Image Processing,
%      12(8):938-943, 2003.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1217270&tag=1>
%
% [ANIGAUSS]  Source code and documentation available at
%      <http://staff.science.uva.nl/~mark/downloads.html>
%
%% See also 
% Related:
% <GAUSSKERNEL.html |GAUSSKERNEL|>,
% <HOURGLASSKERNEL.html |HOURGLASSKERNEL|>,
% <EUCLIDKERNEL.html |EUCLIDKERNEL|>,
% <../../filter/html/CONVOLUTION.html |CONVOLUTION|>.
% Called:
% <matlab:webpub(whichpath('ANIGAUSS')) |ANIGAUSS|>.

%% Function implementation
function H = dirgausskernel(sigma1, sigma2, theta, m, method)

%%
% parsing/checking parameters

error(nargchk(3, 5, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

if ~all(size(sigma1)==size(sigma2)) || ~all(size(sigma1)==size(theta))
    error('dirgausskernel:inputerror', ...
        'input sigmas and theta must have same dimension');
    
elseif strcmp(method,'geu') && ~exist('anigauss','file')
    error('dirgausskernel:unknownfunction','anigauss toolbox needs to be loaded');
end


%%
% main calculation

n = length(sigma1(:));

% size of the kernel
s = max([sigma1(:);sigma2(:)]);
m = min( s*3, m);
m = round(m/2)*2 + 1;

if strcmp(method,'geu'),    
    h = zeros(m,m); 
    c = (m + 1)/2;
    h(c,c) = 1;
end

% create the output filters
H = zeros( m, m, n);   

for i=1:n
    
    s1 = sigma1(i);
    s2 = sigma2(i);
    t = theta(i);
    
    switch method
        
        case 'geu'
            angle =  t / pi *180;
            % trick of ANIGAUSS: note the order of (s1,s2)... only way to
            % respect the orientation
            if s1<s2
                tmp = s1;
                s1 = s2;
                s2 = tmp;
            end
            H(:,:,i) = anigauss(h, s1, s2, angle, 0, 0);
            
        case '2d'
            s1 = s1*s1;
            s2 = s2*s2;
            x = linspace(-m/2,m/2,m);
            [Y,X] = meshgrid(x,x);
            % rotate
            X1 = cos(t)*X + sin(t)*Y;
            Y1 = - sin(t)*X + cos(t)*Y;
            X1 = X1 / s1;
            Y1 = Y1 / s2;
            h = exp( -(X1.^2 + Y1.^2) );
            h = h/sum(h(:));
            H(:,:,i) = h;
            return
        
    end
    
end

end % end of dirgausskernel
