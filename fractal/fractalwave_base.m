function [S, Histoh, H, H_r, DH, ErrDh, varargout] = ...
    fractalwave_base(I, measure, rho, sc0, scw, scr, scn, hmin, hmax, hstep, ...
    sigma, der, int, samp, gn, tn)
% FRACTALWAVE_BASE - Base function for FRACTALWAVE.
% 
%    [S, Histoh, H, H_r, DH, ErrDh] = ...
%      FRACTALWAVE_BASE(I, measure, rho, sc0, scw, scr, scn, ...
%                       hmin, hmax, hstep, sigma, der, int, samp, gn, tn);
%    [S, Histoh, H, H_r, DH, ErrDh, edges, scales, Tpsi, R] = ...
%      FRACTALWAVE_BASE(I, measure, rho, sc0, scw, scr, scn, ...
%                       hmin, hmax, hstep, sigma, der, int, samp, gn, tn);
%
% credit: J.Grazzini (ISR-2/LANL)
% 
% See also 
% Related: FRACTALWAVE, FRACTALMORPH_BASE --
% Called: WAVEPROFILE_BASE, SINGFRACTAL_BASE.

% windex = 0
% ord_der = 0
% sc0 [0][0] = 1
% sw [0][0] = 16
% wavrange = 10
% wavpoints = 7
% scale_wavelet = sw / max(x,y)
% qs = wavrange^(1/(wavpoints-1))
% (sc=1, ip=0; ip<=wavpoints; ip++, sc=sc*qs)
% xx = 1 / log(sc * sc0)
% generate_wavelet_normalize at scale sc
% convolve(T)
% yy = log(abs(T))
% yy *= xx
% log(T) = log(alpha) + h*log(r)

narginchk(7, 16);
nargoutchk(1, 10);

% we allow a variable number of inputs with this function
if nargin<16,  tn = false;
    if nargin<15,  gn = false;
        if nargin<14,  samp = 1;
            if nargin<13,  int = 'fast';
                if nargin<12,  der = 'fast';
                    if nargin<11,  sigma = 0;
                        if nargin<10,  hstep = 100;
                            if nargin<9,  hmax = 2;
                                if nargin<8,  hmin  = -1;  end
                            end
                        end
                    end
                end
            end
        end
    end
end


%% Dealing with multispectral in the case the method is decorrelated

[X,Y,C] = size(I); % possibly multichannel when C>1

if C>1 && any(strcmp(measure,{'grad','scalar','decorr'}))
    S = cell(C,1);
    Histoh = cell(C,1);
    H = cell(C,1);
    H_r = cell(C,1);
    DH = cell(C,1);
    ErrDh = cell(C,1);
    if nargout>=9,  varargout{2} = cell(C,1);
        if nargout==10,  varargout{3} = cell(C,1);  end;
    end
    for c=1:C
        [s, histoh, h, h_r, dh, errDh, Rho, Edges, Tpsi, R] = ...
            fractalwave_base(I(:,:,c), 'grad', rho, ...
            sc0, scw, scr, scn, hmin, hmax, hstep, ...
            sigma, der, int, samp, gn, tn);
        S{c} = s;
        Histoh{c} = histoh;
        H{c} = h;
        H_r{c} = h_r;
        DH{c} = dh;
        ErrDh{c} = errDh;
        if nargout>=9,  varargout{3}{c} = Tpsi;
            if nargout==10,  varargout{4}{c} = R;  end;
        end
    end
    if nargout>=7,  varargout{1} = Rho;  end
    if nargout>=8,  varargout{2} = Edges;  end
    return;     
end

%% Setting internal variables

trans = false;

% WINDEX=0;
% WORD_DER=0;
sc0 = sc0 * scale_wave(scw, X, Y); % output.sc0

if isempty(rho)
    % scale step
    scs = scr ^ (1/(scn-1)); % qs=pow(WAV_RANGE,1./(WPOINTS-1));
    % scales
    rho = scs .^ (0:scn-1); % sc=qs^i, i=0:scn-1

else
    n_rho =  length(rho);
    if n_rho==1
        rho = linspace(sigma, rho, 10);
    elseif n_rho==2
        rho = linspace(rho(1), rho(2), n_rho);
    end
    
end

rho = rho(:); % rho is transformed to a column vector
n_rho =  length(rho);

pad = 1;   %pad=0;

A = I;
%h = ceil(max(6*rr+1,3));
%A = padarray(I,[h h],'replicate','both');

%% Perform wavelet transforms at different scales

Tpsi = waveprofile_base(A, rho, measure, sigma, der, int, samp, gn, tn);

if pad  
    nx = X-2*pad; ny = Y-2*pad;                                        
    Tpsi = Tpsi(pad+1:end-pad, pad+1:end-pad, :);
else
    nx = X; ny = Y;                                                    %#ok
end
n_x = nx * ny;

% transform the matrix by permuting/reshaping its column
Tpsi_r = permute(reshape(Tpsi, [n_x n_rho]), [2 1]);
% Tpsi_r is a 2D matrix which stores on the i-th column all the wavelet
% coefficients computed on n_rho scales for the i-th pixel


%% Compute the singularity exponents as the regression coefficients

[S, R, Histoh, H, H_r, DH, ErrDh, edges] = ...
    singfractal_base(Tpsi_r, rho, sc0, hmin, hmax, hstep, trans);

% reshape the matrix of singularity exponents and the matrix of correlation
% coefficients
S = reshape(S, [nx ny]);
R = reshape(R, [nx ny]);

%% Return the desired outputs

if nargout>=7,  varargout{1} = [rho; sc0];
    if nargout>=8,  varargout{2} = edges;
        if nargout>=9,  varargout{3} = Tpsi;
            if nargout>=10,  varargout{4} = R;
            end
        end
    end
end

end


%--------------------------------------------------------------------------
function sc0 = scale_wave(scw, X, Y)
sc0 = scw / max(X,Y);
% sc0 = scw /  sqrt(X*Y);
end



