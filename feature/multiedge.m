function map = multiedge(I,nsc,varargin)
% MULTIEDGE - Compute a multiscale edge map using a 2D Sobel-like wavelet.
%
%         map = multiedge(I,nsc);
%         map = multiedge(I,nsc,ave);
%
% Inputs:
%     I : original image.
%     nsc : maximum level of wavelet scales.
%     ave : size of the filter used
%
% Output:
%     map : binary edge map.

% acknowledgment: "A Wavelet Approach to Edge Detection" by J.Li, 2003
%    http://eref.uqu.edu.sa/files/Others/Image%20Processing/A%20wavelet%20a
%    pproach%20to%20edge%20detection%20-%20Thesis.pdf

%% Parsing parameters

p = inputParser;   % create an instance of the inputParser class.
% mandatory parameter
p.addRequired('I', @isnumeric);
p.addRequired('nsc', @(x)isscalar(x) && round(x)==x && x>=1);
p.addOptional('ave', 10, @(x)isscalar(x) && round(x)==x && x>1);

% parse and validate all input arguments
p.FunctionName='MULTIEDGE';
p.parse(I,nsc,varargin{:}); 

% create the variables associated to the fieldnames
ave = p.Results.ave;

%% Main computation

[m n] = size(I);

h1 = fspecial('sobel');
h2 = h1';
v = [1 2 1;2 4 2; 1 2 1];

I1 = filter2(h1,I);
I2 = filter2(h2,I);

ns = h1;
% loop over the scales range
for ii=1:nsc
    nsh = conv2(ns,v); % imfilter(ns,v);
    p = max(max(ns)) / max(max(nsh));
    
    % increase the scale of the wavelet
    Ih1 = filter2(v,I1) * p;
    Ih2 = filter2(v,I2) * p;
    ns = nsh;
    
    % distinguish noise and real edges
    I1 = I1.*(abs(I1)<=abs(Ih1)) + Ih1.*(abs(I1)>abs(Ih1));
    I2 = I2.*(abs(I2)<=abs(Ih2)) + Ih2.*(abs(I2)>abs(Ih2));
end

% compute the magnitude
e = (I1.^2 + I2.^2);

% threshold the magnitude
ave = fspecial('average', ave);

map = (e>(imfilter(e,ave))) .* (e>mean2(e));
end


%--------------------------------------------------------------------------
function y = psi(axis,a,h) 
% PSI - cascade algorithm for the wavelet model 
%
%        y = psi(axis, a, h);
%
% Inputs:
%    axis : 1 for x, 2 for y, 3 for (x+y), 4 for (x-y) 
%    a : variance 
%    h : order of derivative 
% Output:
%    y : wavelet with support = (a+h+1)x(a+h+1) 
%
% Examples
%   y=psi(1,10,0) gives a Gaussian with support of 11 by 11 
%   y=psi(3,6,4) gives a 4th derivative of Gaussian rotated pi/4. 

switch axis
    case 'x'
        d = [0 0 0; 0 1 -1; 0 0 0]/2;
    case 'y'
        d = [0 0 0; 0 1 0; 0 -1 0]/2;
    case {'x+y','y+x'}
        d = [0 0 0; 0 1 0; 0 0 -1]/2;
    case 'x-y'
        d = [0 0 0; 0 1 0; -1 0 0]/2;
    otherwise
        error('psi:invalidparameter','invalid axis parameter');
end

v = [1 1; 1 1]/4;

if h > 0
    y = d;
    for i=1:h-1,  
        y = conv2(y,d);
    end
    
else
    y = v;
    a = a-1;
end

for i=1:a
    y = conv2(y,v);
end

end