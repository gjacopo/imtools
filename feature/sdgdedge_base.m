%% SDGDEDGE_BASE - Base function for SDGDEDGE.
%
%% Syntax
%     edgemap = SDGDEDGE_BASE(I, plus, sigma, mu, thres, hsize);
%
%% Acknowledgment:
% <mailto:sergei.koptenko{at}resonantmedical.com  Sergei Koptenko>, Resonant
% Medical, Montreal (Qc., Canada), <www.resonantmedical.com>.  
%           
%% See also 
% Related:
% <SDGDEDGE.html |SDGDEDGE|>,
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>,
% <CANNYEDGE_BASE.html |CANNYEDGE_BASE|>,
% <CANNYEDGEMAP_BASE.html |CANNYEDGEMAP_BASE|>,
% <ROTHWELLEDGE_BASE.html |ROTHWELLEDGE_BASE|>,
% <CONGRUENCYEDGE_BASE.html |CONGRUENCYEDGE_BASE|>,
% <COMPASSEDGE_BASE.html |COMPASSEDGE_BASE|>,
% <ANISOEDGE_BASE.html |ANISOEDGE_BASE|>,
% <ELDERZUCKEREDGE_BASE.html |ELDERZUCKEREDGE_BASE|>,
% <KOETHEDGE_BASE.html |KOETHEDGE_BASE|>,
% <PETROUEDGE_BASE.html |PETROUEDGE_BASE|>.
% Called:
% <matlab:web(whichpath('CONV2')) |CONV2|>.

%% Function implementation
function edgemap = sdgdedge_base(I, plus, sigma, mu, thres, hsize)

%% 
% internal parameters
C = size(I,3);
if isempty(hsize),    hsize = ceil(3 * sigma);  end;
x = 2*hsize + 1;


%% 
% dealing with multispectral images
if C>1
    edgemap = false(size(I));
    for c=1:C
        edgemap(:,:,c) = sdgdedge_base(I(:,:,c), plus, sigma, mu, thres, hsize);
     end
    return;
end

%%
% compute the derivarives kernels
dGx = GaussDx2(x, mu, sigma, hsize); % first derivative dGx
dGy = dGx';                             % first derivative dGy
dGxx = conv2(dGx, dGx, 'same');         % second Derivative dGxx
dGyy = dGxx';                           % second Derivative dGyy
dGxy = conv2(dGx, dGy, 'same');         % second Derivative dGxy 

%% 
% estimate the derivatives
a_x = conv2(I, dGx, 'same');
a_y = conv2(I, dGy, 'same'); 
a_xx = conv2(I, dGxx, 'same');
a_yy = conv2(I, dGyy, 'same');
a_xy = conv2(I, dGxy, 'same');
% warning off MATLAB:divideByZero 
% Warning: derivatives a_x and a_y for some pixels may be zero; this produce 
% an error "divideByZero" and such pixels gets value =NaN. It is recomended
% to check for NaNs after filtering and replacing it with any value of 
% convenience: zero, global minimum, etc...  
% warning on MATLAB:divideByZero

%%
% output of the filter
edgemap = (a_xx .* a_x.^2 + 2* a_xy .* a_x .* a_y + a_yy .* a_y .^2 );
edgemap = edgemap ./ (a_x .^2 + a_y .^2);

if plus,      edgemap =  edgemap + a_xxm + a_yy;   end

%%
% final binary map
mm = max(edgemap(:));
edgemap = edgemap / mm > thres;

end % end of sdgdedge_base


%% Subfunction

%%
% |GAUSSDX2| - Compute a 2D Gaussian Derivative.
%
% Inputs:
%   |gdsize| : size of Gaussian kernel.
%
%   |mu| : mean of the Gaussian.
%
%   |sigma| : standard deviation of a Gaussian function.
%
%   |sigma_width| : defines where to cut the Gaussian kernel tail (or width 
%      of the kernel in sigma); default: |sigma_width=3| or 98% of Gaussian.
%
% Output:
%   |gausskernel_2d| : derivative kernel.
%--------------------------------------------------------------------------
function GaussKernel_2D = GaussDx2(GDsize, mu, sigma, sigma_width)

x = linspace(-sigma_width, sigma_width, GDsize); 

Gaussian = exp(- 0.5*(((x - mu)/ sigma) .^2)) /(sigma * realsqrt(2 * pi)); 
Gaussian = Gaussian / sum(Gaussian); % normalize the sum of samples to 1

GaussKernel_1D = (mu - x) .* Gaussian / sigma^2;
% normalize kernel
GaussKernel_1D  = GaussKernel_1D  /sum(GaussKernel_1D .* (mu-x)); 
GaussKernel_2D = conv2(GaussKernel_1D, Gaussian', 'full');
end % end of GaussDx2
