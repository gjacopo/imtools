function [S, varargout] = fractalwave(I, varargin)
% FRACTALWAVE.
% 
%      S = FRACTALWAVE(I);
%      S = FRACTALWAVE(I, measure, rho);
%      [S, Histoh, H, H_r, DH, ErrDh, edges, scales, Tpsi, R ] = ... 
%         FRACTALWAVE(I, measure, rho, 'Property', propertyvalue, ... );
%
% Inputs:
%   I :
%   measure :  default: =; 
%   rho :  default: rho=[];
% 
% Property [propertyname  propertyvalues]
%   'sc0' :  default: sc0=1;
%   'scw' :  default: scw=16;
%   'scr' :  default: scr=10;
%   'scn' :  default: scn=7;
%   'hmin' :  default: hmin=-1;
%   'hmax' :  default: hmax=2;
%   'hstep' :  default: hstep=100;
%   'hcorr' :  default: hcorr=0.999;
%   'sig' :  default: sig=0;
%   'der' :  default: der='fast';
%   'int' :  default: int='fast';
%   'samp' :  default: samp=1;
%   'tn' :  default: tn=false;
%   'gn' :  default: gn=false;
%
% Outputs:
%   S :
%   Histoh :
%   H :
%   H_r :
%   DH :
%   ErrDh :
%   scales :
%   edges :
%   Tpsi :
%   R :
%   
% credit: J.Grazzini (ISR-2/LANL)
%
% See also 
% Related: FRACTALWAVESTAT, FRACTALMORPH --
% Called: FRACTALWAVE_BASE.

%% Parsing parameters

narginchk(1, 37);
nargoutchk(0, 10);

% mandatory parameter
if ~isnumeric(I)
    error('fractalwave:inputerror','a matrix is required in input'); 
end

% optional parameters
p = createParser('FRACTALWAVE');   

% additional parameters regarding numerical estimation
p.addOptional('measure', 'tens', @(x)ischar(x) && ...
    any(strcmpi(x,{'grad','scalar','decorr','tens','riemann','corr'})));

% principal scale parameters
p.addOptional('rho', [], @(x)isfloat(x) && length(x)>=1 && any(x)>=0);

% options for multiscale analysis
p.addParamValue('sc0', 1, @(x)isscalar(x) && x>0);  % output.sc0
p.addParamValue('scw', 16, @(x)isscalar(x) && x>0); % SC_W[0][0]
p.addParamValue('scr', 10, @(x)isscalar(x) && x>0); % WAV_RANGE
p.addParamValue('scn', 7, @(x)isscalar(x) && x>1);  % WPOINTS
% minimum and maximum singularity values considered in the histogram
p.addParamValue('hmin', -1, @(x)isscalar(x) && isfloat(x) && x<=0);
p.addParamValue('hmax', 2, @(x)isscalar(x) && isfloat(x) && x>=0 );
p.addParamValue('hstep', 100, @(x)isscalar(x) && ...
    (isfloat(x) && x>0 && x<=1) || (x>1 && round(x)==x));
% goodness for singularity analysis regression
p.addParamValue('hcorr', 0.999, @(x)isscalar(x) && isfloat(x) && x>=0 && x<=1);
% options for estimating the GST
p.addParamValue('sig', 0, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addParamValue('der', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','vista','fast','conv','fleck', ...
    'tap5','tap7','sob','opt','ana'}))));
p.addParamValue('int', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','conv','fast','ani'}))));
p.addParamValue('samp', 1, @(x)isscalar(x) && round(x)==x && x>=1 && x<=5);
p.addParamValue('gn', false, @(x)islogical(x));
p.addParamValue('tn', false, @(x)islogical(x));

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

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            


%% Checking/setting parameters

[~,~,C] = size(I);
%n_x = X*Y;


%% Compute the multifractal decomposition

[S, Histoh, H, H_r, DH, ErrDh, scales, edges, Tpsi, R] = ...
    fractalwave_base(I, p.measure, p.rho, p.sc0, p.scw, p.scr, p.scn, ...
    p.hmin, p.hmax, p.hstep, p.sig, p.der, p.int, p.samp, p.gn, p.tn);
% note: the last element of scales is sc0

if nargout>=2, varargout{1} = Histoh;
    if nargout>=3, varargout{2} = H;
        if nargout>=4, varargout{3} = H_r;
            if nargout>=5, varargout{4} = DH;
                if nargout>=6,  varargout{5} = ErrDh;
                    if nargout>=7, varargout{6} = scales;
                        if nargout>=8, varargout{7} = edges;
                            if nargout>=9,  varargout{8} = Tpsi;
                                if nargout==10,  varargout{9} = R;  end
                            end
                        end
                    end
                end
            end
        end
    end
end


%% Display

if p.disp
        
    ndisp = 3; mdisp = ceil((length(scales)-1)/ndisp);
    
    % we distinguish the cases 'riemann' and 'scalar'
    if C>1 && any(strcmp(p.measure,{'grad','scalar','decorr'}))
        n_x = numel(S{1});
        for c=1:C
            % display the wavelet transform
            figure % open a display window
            for r=1:length(scales)-1
                subplot(mdisp, ndisp, r),  imagesc(rescale(Tpsi{c}(:,:,r))),
                colormap gray, title(['WT at \sigma=' num2str(scales(r))])
            end
            S{c}(S{c}>p.hmax) = p.hmax;
            h0 = min(S{c}(:));     h1 = max(S{c}(:));
            %delta_h = 0.01;      range_h = (h1 - h0) / delta_h;
            % output range of singularities
            disp_sing = ['Singularity range for channel ' num2str(c)];
            disp([disp_sing ': [ ' num2str(h0) ' - ' num2str(h1) ' ]']);
            figure, % open another window
            % display the singularity exponents
            title_sing = ['singularity exponents for channel ' num2str(c)];
            subplot(1,2,1), imagesc(S{c}), axis image off, colormap gray,
            title(title_sing);
            % output regression rate
            disp_reg = ['Good regression points for channel ' num2str(c)];
            disp([disp_reg ': ' num2str(100*sum(abs(R{c}(:))>=p.hcorr)/n_x) ' %']);
            subplot(1,2,2),imagesc(abs(R{c})>=p.hcorr), axis image off, colormap gray,
            title('''good regression'' map');
            % display its histogram
            title_exp = ['exponents distribution for channel ' num2str(c)];      
            figure, bar(edges,Histoh{c},'histc'), % hist(S{c}(:),range_h)
            title(title_exp)
            title_spec = ['multifractal spectrum for channel ' num2str(c)];
            figure, errorbar(H_r{c}, DH{c}, ErrDh{c}/2, ErrDh{c}/2), title(title_spec);
        end
        
    else
        n_x = numel(S);
        figure;
        for r=1:length(scales)-1
            subplot(mdisp, ndisp, r),  imagesc(rescale(Tpsi(:,:,r))),
            colormap gray, title(['WT at \sigma=' num2str(scales(r))]);
        end
        S(S>p.hmax) = p.hmax;
        h0 = min(S(:));    h1 = max(S(:));
        disp(['Singularity range: [ ' num2str(h0) ' - ' num2str(h1) ' ]']);
        figure, subplot(1,2,1), imagesc(S), axis image off, colormap gray, 
        title('singularity exponents');
        disp(['Good regression points: ' ...
            num2str(100*sum(abs(R(:))>=p.hcorr) / n_x) ' %']);
        subplot(1,2,2), imagesc(abs(R)>=p.hcorr), axis image off, colormap gray,
        title('''good regression'' map');
        figure, bar(edges,Histoh,'histc'), title('exponents distribution');
        figure, errorbar(H_r, DH, ErrDh/2, ErrDh/2), title('multifractal spectrum');
    end
end

end



