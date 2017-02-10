function [S, R, Histoh, H, H_r, DH, ErrDh, edges] = ...
    singfractal_base(MS, scale, sc0, hmin, hmax, hstep, trans)
% SINGFRACTAL_BASE - Compute the multifractal singularity exponents and the
% multifractal spectrum.
%
%       [S, R, Histoh, H, H_r, DH, ErrDh, edges] = ...
%           SINGFRACTAL_BASE(MS, scale, sc0, hmin, hmax, hstep, trans)
%
% credit: J.Grazzini (ISR-2/LANL)
%
% See also 
% Related: FRACTALWAVE_BASE, SINGFRACTALRECONS_BASE, FRACTALMORPH_BASE.

narginchk(3, 7);
nargoutchk(1, 8);

% we allow a variable number of inputs with this function
if nargin<7,  trans = false;
    if nargin<6,  hstep = 100;
        if nargin<5,  hmax = 2;
            if nargin<4,  hmin  = -1;  end
        end
    end
end

[S,R] = calculate_multifractal(MS, scale, sc0);

[Histoh, H, edges] = singularity_histogram(S, hmin, hmax, hstep, sc0, trans);

[H_r, DH, ErrDh] = spectrum_histogram(H, Histoh, sc0);

end



%--------------------------------------------------------------------------
function [S,R] = calculate_multifractal(Tpsi_r, rho, sc0)

funreg = @(x)log(x); 

n_x = size(Tpsi_r,2);

% compute the average
Tpsi_mean = mean(Tpsi_r);                                              %#ok

% scale variable considered for the regression: log(r)
rr =  1 ./ funreg(rho*sc0); % xx[ip] = 1 / log(sc*output.sc0);
% rr is transformed for estimating the log of the wavelet transform divided
% by the log of the scale
rr = repmat(rr,[1 n_x]); % the entry rr is repeated over the columns
% compute the transform log(|T_psi|)
Tpsi_r = funreg(abs(Tpsi_r));   % yy[ip] = log(fabs(wave[ip]));
% compute the transform log(|T_psi|)/log(r): yy[ip] *= xx[ip];
Tpsi_r = Tpsi_r .* rr;

% compute the regression coefficients=singularity exponents for each pixel
[~, S, R] = linearfit(rr, Tpsi_r);
% returns, for each pixel, the singularity exponent in the intercept b of
% the log-log linear regression

end


%--------------------------------------------------------------------------
function [a, b, R, E] = linearfit(xx, yy)
% linear regression - fast vectorized implementation

sumx = mean(xx);
sumy = mean(yy);

sx = mean(xx.^2);
sxy = mean(xx.*yy);
sy = mean(yy.^2);

sxy = sxy - sumx.*sumy;
sx = sqrt(abs(sx - sumx.*sumx));
sy = sqrt(abs(sy - sumy.*sumy));

I=sx>eps;

a(I) = sxy(I) ./ (sx(I).*sx(I));
b(I) = sumy(I) - a(I) .* sumx(I);
J=sy>eps;
R(I&J) = sxy(I&J) ./ (sx(I&J) .*sy(I&J));
R(I&~J) = 1;

a(~I) = 0;
b(~I) = sumy(~I);
R(~I) = 1;

E = sy .* sy - a .* sx .*sx;

% non vectorized matlab implementation using regress, polyfit or \
% n_x = size(yy,2);
% for x=1:n_x
%     % \ : Backslash or left matrix divide.
%     % A\B is the matrix division of A into B. If A is an N-by-N matrix and
%     % B is a matrix with several such columns, then X = A\B is the solution
%     % to the equation A*X = B.
%     xx = xx(:,1); % xx is 'back'-transformed for using the operator \
%     %XX = [ones(length(xx),1), xx];
%     rr = xx \ yy(:,x);
%     a(x) = rr(2);  b(x) = rr(1);                                 
%     % p = polyfit(xx,yy(:,x),1);
%     % a(x) = p(1); b(x) = p(2);   
%     % p1 = regress(yy,xx);
%     % a(x) = p1(2); b(x) = p1(1);  
%     [rr,P] = corrcoef([xx yy(:,x)]);
%     R(x) = rr(1,2);
% end

end


%--------------------------------------------------------------------------
function [histoh, h, edges] = ...
    singularity_histogram(S, hmin, hmax, hstep, sc0, trans)

S = S(:);

%h0 = min(S);    h1 = max(S);                                           
%delta_h = (h1 - h0) / hstep;

if hstep>0 && hstep<1
    delta_h = hstep;
elseif hstep>1
    delta_h = (hmax - hmin) / hstep;                                   
end
    
S(S>hmax | isnan(S)) = hmax+0.1; % S(S>hmax) = [];
S(S<hmin) = hmin-0.1;

% histoh = hist(S(:),range_h);
% [histoh,bin] = histc(S(:),h0:delta_h:h1);
edges = [hmin-0.1 hmin:delta_h:hmax hmax+0.1];
[histoh,bins] = histc(S, edges);
% note: with the previous transformation of S, there should not be 'out of
% range' values (ie. null values in bins) anymore
h = accumarray(bins, S, [length(histoh) 1]);

% get rid of those pixels outside the bounds [hmin,hmax]
histoh = histoh(2:end-1);
h = h(2:end-1);
% closing for the max: count/store the values for which S=hmax with the bin
% hmax-delta_h<=S<hmax
histoh(end-1) = histoh(end-1) + histoh(end); histoh(end) = [];
h(end-1) = h(end-1) + h(end); h(end) = [];
% +outliers

% also modify edges accordingly
edges = edges(2:end-2);

if trans
    [hmode, ihmode] = max(histoh);
    
    I = histoh(1:ihmode)>0;
    auxdist = h(I) ./ histoh(I) + log(histoh(I)/hmode) / log(sc0);
    hshift = min(auxdist);
    
    h = h - hshift * histoh;
end

end


%--------------------------------------------------------------------------
function [h_r, Dh, errDh] = spectrum_histogram(h, histoh, sc0)

min_ev = 100;
cump = cumsum(histoh);
%I = floor(cump/min_ev);
%[dum,ir0] = unique(I,'first');
%[dum,ir1] = unique(I,'last');

% I = [0; diff(floor(cump/min_ev))]
% ir0 = find(I);
% ir1 = find(~I);
% if length(ir1)<length(ir0),  ir1 = [ir1; ir0(end)];  end

Nh = length(histoh);
histo_r = zeros(Nh,1);
h_r = zeros(Nh,1);

i = 1; ir0 = 1;
ir1 = find(cump>=min_ev,1,'first');
while ir0<=length(histoh) && ~isempty(ir1)
    histo_r(i) = sum(histoh(ir0:ir1).^2) / cump(ir1); % sum(histoh(ir0:ir1));
    h_r(i) = sum(h(ir0:ir1)) / cump(ir1);
    cump = cump - cump(ir1); cump(1:ir1) = 0;
    ir0 = ir1 + 1;
    ir1 = find(cump>=min_ev,1,'first');
    i = i+1;
end
if ~all(cump==0)
    histo_r(i) = sum(histoh(ir0:end).^2) / cump(end); % sum(histoh(ir0:ir1));
    h_r(i) = sum(h(ir0:end)) / cump(end);
end

histo_r(i+1:end) = []; 
h_r(i+1:end) = [];

% the confidence range is taken as +-Ks sigmas in the distribution of
% probability boxes, which is renormalized binomial. we directly propagate
% the confidence range to the D(h) calculation
Ks = 3;
prob = Ks * sqrt((1-histo_r/sum(histo_r)) ./ histo_r);
prob(prob<eps) = 0;
prob(prob>=1) = 1;

% the error bar is constructed by propagation: as the propagated interval
% is symmetric, we take the bar as the maximum of the two distances to the
% central value. this is given by the lower bound in the logarithm
errDh = log(1-prob) / log(sc0);
errDh(isinf(errDh)) = 1; % cases where prob==1

% finding and normalizing by the mode
maxprob = max(histo_r);
if maxprob>eps,  histo_r = histo_r / maxprob;  end;

% evaluating experimental reduced D(h): D_reduced(h) = D(h) - D_space
histo_r(histo_r==0) = eps;
Dh = - log(histo_r) / log(sc0);

end
