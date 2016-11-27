%% COMPASSEDGE_BASE - Base function for COMPASSEDGE.
%
%% Syntax
%     [S, O, A, U] = COMPASSEDGE_BASE(I, R, Angles, wedges, gray);
%
%% Acknowledgment
% This function calls the mex functions implemented in the original [COMPASS]
% code by Ruzon & Tomasi (available at <http://ai.stanford.edu/~ruzon/compass/>).
%
%% See also
% Related:
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>,
% <CANNYEDGE_BASE.html |CANNYEDGE_BASE|>,
% <ROTHWELLEDGE_BASE.html |ROTHWELLEDGE_BASE|>,
% <CONGRUENCYEDGE_BASE.html |CONGRUENCYEDGE_BASE|>,
% <ANISOEDGE_BASE.html |ANISOEDGE_BASE|>,
% <ELDERZUCKEREDGE_BASE.html |ELDERZUCKEREDGE_BASE|>,
% <KOETHEDGE_BASE.html |KOETHEDGE_BASE|>,
% <SDGDEDGE_BASE.html |SDGDEDGE_BASE|>,
% <PETROUEDGE_BASE.html |PETROUEDGE_BASE|>.
% Called: 
% GREYCOMPASS_MEX, COMPASS_MEX.

%% Function implementation
function [S, O, A, U] = compassedge_base(I, R, Angles, wedges, gray)

%%
% dealing with multichannel images

C = size(I,3);

if C>1 && (C~=3 || (C==3 && islogical(gray) && gray))
    % process channel by channel
    S = zeros(size(I));
    if nargout>=2, varargout{1}= cell(C,1);
        if nargout>=3, varargout{2} = zeros(size(I));
            if nargout==4, varargout{3} = cell(C,1); end
        end
    end
    for c=1:C
        [S(:,:,c), O, A, U] = ...
            compassedge_base(I(:,:,c), R, Angles, wedges, gray);
        if nargout>=2, varargout{1}{c} = O;
            if nargout>=3, varargout{2}(:,:,c) = A;
                if nargout==4, varargout{3}{c} = U; end
            end
        end
    end  
    return;
    
end

%%
% run the mex files

maxradius = ceil(3 * R); % see functions compass_mex.c & compass.c
A = padarray(I,[maxradius maxradius],'replicate','both');

% check  again the size of the input image
C = size(I,3);

% default parameters
spc = 1; % spacing at each scale
Is = 0; % the compass operator is applied on the entire image.

if C==3
    [S, O, A, U] = compass_mex(A, R, spc, Is, Angles, wedges);

else % note that we have necessarly C=1 at this stage
     [S, O, A, U] = greycompass_mex(A, R, spc, Is, Angles, wedges);
end

%%
% refine outputs

S = S(1:end-1,1:end-1);
A = A(1:end-1,1:end-1);
O = O(1:end-1,1:end-1,:);
U = U(1:end-1,1:end-1,:);

end % end of compassedge_base
