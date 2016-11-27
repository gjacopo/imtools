%% CONGRUENCYEDGE_BASE - Base function for CONGRUENCYEDGE.
%
%% Syntax
%     M  = CONGRUENCYEDGE_BASE(I, nscale, norient, minWaveLength, ...
%                mult, sigmaOnf, k, cutOff, g, noiseMethod);
%     [M m or ft pc EO T]  = CONGRUENCYEDGE_BASE(I, nscale, norient, ...
%                minWaveLength, mult, sigmaOnf, k, cutOff, g, noiseMethod);
% 
%% Acknowledgement
% This is just a wrapper to P.Kovesi's edge/corner detectors based on phase
% congruency. See the webpage:
%       <http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/index.html>
% In particular, it enables to deal with multispectral images.
%
%% See also
% Related:
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>,
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>,
% <CANNYEDGE_BASE.html |CANNYEDGE_BASE|>,
% <ROTHWELLEDGE_BASE.html |ROTHWELLEDGE_BASE|>,
% <COMPASSEDGE_BASE.html |COMPASSEDGE_BASE|>,
% <ANISOEDGE_BASE.html |ANISOEDGE_BASE|>,
% <ELDERZUCKEREDGE_BASE.html |ELDERZUCKEREDGE_BASE|>,
% <KOETHEDGE_BASE.html |KOETHEDGE_BASE|>,
% <SDGDEDGE_BASE.html |SDGDEDGE_BASE|>,
% <PETROUEDGE_BASE.html |PETROUEDGE_BASE|>.
% Called: 
% <../../../toolboxes/kovesi/PhaseCongruency/html/PHASECONG3.html |PHASECONG3|>.
 
%% Function implementation
function [M, varargout]  = ...
    congruencyedge_base(I, nscale, norient, minWaveLength, ...
    mult, sigmaOnf, k, cutOff, g, noiseMethod)

%%
% dealing with multispectral images

C = size(I,3);

if C>1
    M = zeros(size(I));
    if nargout>=2,
        varargout{1} = zeros(size(I)); % m
        if nargout>=3,
            varargout{2} = zeros(size(I)); % or
            if nargout>=4,
                varargout{3} = zeros(size(I)); % ft
                if nargout>=5,
                    varargout{4} = cell(C,1); % pc
                    if nargout>=6,
                        varargout{5} = cell(C,1); % EO
                        if nargout==7,
                            varargout{6} = zeros(size(I)); % T
                        end
                    end
                end
            end
        end
    end
    for c=1:C
        [M(:,:,c) m or ft pc EO, T]  = ...
            congruencyedge_base(I(:,:,c), nscale, norient, minWaveLength, ...
            mult, sigmaOnf, k, cutOff, g, noiseMethod);
        if nargout>=2,
            varargout{1}(:,:,c) = m;
            if nargout>=3,  
                varargout{2}(:,:,c) = or;
                if nargout>=4,  
                    varargout{3}(:,:,c) = ft;
                    if nargout>=5,  
                        varargout{4}{c} = pc;
                        if nargout>=6,  
                            varargout{5}{c} = EO;
                            if nargout==7,  
                                varargout{6}(:,:,c) = T;
                            end
                        end
                    end
                end
            end
        end
    end
    return
end

%%
% main calculation: call to |PHASECONG3|

[M m or ft pc EO, T] = phasecong3(I, nscale, norient, minWaveLength, ...
    mult, sigmaOnf, k, cutOff, g, noiseMethod);

if nargout>=2,
    varargout{1} = m;
    if nargout>=3,
        varargout{2} = or;
        if nargout>=4,
            varargout{3} = ft;
            if nargout>=5,
                varargout{4} = pc;
                if nargout>=6,
                    varargout{5} = EO;
                    if nargout==7,
                        varargout{6} = T;
                    end
                end
            end
        end
    end
end

end % end of congruencyedge_base

