function [Histoh, H, H_r, DH, errDH] = ...
    fractalwavestat(coll, measure, rho, sc0, scw, scr, scn, hmin, hmax, hstep, ...
    sigma, der, int, samp, gn, tn)
% FRACTALWAVESTAT
%
%    Histoh, H, H_r, DH, errDH] = ...
%        FRACTALWAVESTAT(coll, measure, rho, sc0, scw, scr, scn, ...
%                        hmin, hmax, hstep, sigma, der, int, samp, gn, tn);    
%
% credit: J.Grazzini (ISR-2/LANL)
% 
% See also 
% Related: FRACTALWAVE --
% Called: FRACTALWAVE_BASE.

% we allow a variable number of inputs with this function
if nargin<16,  tn = false;
    if nargin<15,  gn = false;
        if nargin<14,  samp = 1;
            if nargin<13,  int = 'fast';
                if nargin<12,  der = 'fast';
                    if nargin<11,  sigma = 0;
                        if nargin<10,  hstep = 0.05;
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

ncoll = size(coll,1);

for i=1:ncoll
    disp(['processing image #' num2str(i) ': ' strtrim(coll(i,:)) '...']);
    % read the input image
    I = double(imread(strtrim(coll(i,:))));
   
    % process the loaded image
    disp('          computing multifracal representation...');
    [s, histoh, h, h_r, Dh, errDh, scales, edges] = ...
        fractalwave_base(I, measure, rho, sc0, scw, scr, scn, ...
        hmin, hmax, hstep, sigma, der, int, samp, gn, tn);
    
    if i==1, SC0 = scales(end);  end  % initialize with the first loaded image
    
    if i>1,  disp('          merging statistics...');  end
    
    [histoh, h] = histoh_adapt(histoh, h, SC0, scales(end));
    
    if i==1  % initialize with the first loaded image
        Histoh = histoh;
        H = h;
    else
        Histoh = Histoh + histoh;
        H = H + h;
    end
end

% figure, bar(edges,H);
figure, bar(edges,Histoh);


end


%--------------------------------------------------------------------------
function [histoh, h, uh, ah] = histoh_adapt(histoh, h, sc0, sc1)

[hmode, ihmode] = max(histoh);                                         %#ok

% we modify the histogram taking into account the fact that
%    histo(h) = hmode * r^{d-D(h)}
% and now r must change fomr sc1 (that of data) to sc0 (new basic
% resolution). The average h associated to the bin box must also be changed
% accordingly
I = histoh>0.5;
Nnew = hmode * exp(log(histoh(I)/hmode) * log(sc0) / log(sc1));

% even if this rule for scale transformation implicitely use the existence
% of the singularity spectrum, this rule is more general and is valid for
% any quantity exhibiting a general power-law scaling in the histogram. In
% addition, notice that if sc0=sc1, then Nnew=histoh(ih).
h(I) = (h(I) .* Nnew) ./ histoh(I);
histoh(I) = Nnew;

% uniform representative values of h
uh = hmin + delta_h * ((1:length(h))' - 0.5);
% average representative values of h
ah = h ./ histoh;
I = isnan(ah);  ah(I) = uh(I);

end



