%% FRACTALMORPH_BASE - Base function for FRACTALMORPH.
%
%% Syntax
%     [S, Histoh, H, H_r, DH, ErrDh, varargout] = ...
%            FRACTALMORPH_BASE(I, measure, rho, sc0, scw, scr, scn, ...
%                              hmin, hmax, hstep, shape);
%
%% See also
% Related:
% <FRACTALMORPH.html |FRACTALMORPH|>,
% <../../fractal/mftensor/html/FRACTALWAVE_BASE.html |FRACTALWAVE_BASE|>.
% Called:
% <FLATSTREL.html |FLATSTREL|>,
% <GRANULOMETRY_BASE.html |GRANULOMETRY_BASE|>,
% <../../fractal/mftensor/html/SINGFRACTAL_BASE.html |SINGFRACTAL_BASE|>.

%% Function implementation
function [S, Histoh, H, H_r, DH, ErrDh, varargout] = ...
    fractalmorph_base(I, measure, rho, sc0, scr, scn, hmin, hmax, hstep, shape)

error(nargchk(7, 11, nargin, 'struct'));
error(nargoutchk(1, 10, nargout, 'struct'));

% we allow a variable number of inputs with this function
if nargin<10,  hstep = 100;
    if nargin<9,  hmax = 2;
        if nargin<8,  hmin  = -1;  end
    end
end

%%
% dealing with multispectral in the case the method is decorrelated

[X,Y,C] = size(I); % possibly multichannel when C>1

if C>1
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
            fractalmorph_base(I(:,:,c), measure, rho, ...
            sc0, scw, scr, scn, hmin, hmax, hstep, shape);
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

%% 
% setting internal variables

trans = false;

% WINDEX=0;
sc0 = sc0 / max(X,Y);

if isempty(rho)
    % scale step
    scs = scr ^ (1/(scn-1)); % qs=pow(WAV_RANGE,1./(WPOINTS-1));
    % scales
    rho = ceil(scs .^ (0:scn-1)); % sc=qs^i, i=0:scn-1

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

A = I;
%h = ceil(max(6*rr+1,3));
%A = padarray(I,[h h],'replicate','both');


%% 
% perform morphological profile transforms at different scales

se = cell(n_rho,1);
for i=1:n_rho
    se{i} = flatstrel(shape, ceil(rho(i)));
end

Tpsi = granulometry_base( A, measure, se );

Tpsi = reshape(cell2mat(Tpsi), [X Y n_rho]); % we wouldn't need to reshape it if it was not for pad

pad = 1;   %pad=0;
if pad  
    nx = X-2*pad; ny = Y-2*pad;                                        
    Tpsi = Tpsi(pad+1:end-pad, pad+1:end-pad, :);
else
    nx = X; ny = Y;                                                    %#ok                                                
end
n_x = nx * ny;

% transform the matrix by permuting/reshaping its column
Tpsi_r = permute(reshape(Tpsi, [n_x n_rho]), [2 1]);
% Tpsi_r is a 2D matrix which stores on the i-th column all the transform
% coefficients computed on n_rho scales for the i-th pixel


%%
% compute the singularity exponents as the regression coefficients

[S, R, Histoh, H, H_r, DH, ErrDh, edges] = ...
    singfractal_base(Tpsi_r, rho, sc0, hmin, hmax, hstep, trans);

% reshape the matrix of singularity exponents and the matrix of correlation
% coefficients
S = reshape(S, [nx ny]);
R = reshape(R, [nx ny]);

%% 
% return the desired outputs

if nargout>=7,  varargout{1} = [rho; sc0];
    if nargout>=8,  varargout{2} = edges;
        if nargout>=9,  varargout{3} = Tpsi;
            if nargout>=10,  varargout{4} = R;
            end
        end
    end
end

end % end of fractalmorph_base
