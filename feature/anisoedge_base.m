%% ANISOEDGE_BASE - Base function for ANISOEDGE.
%
%% Syntax
%     [edgemap, mag] = ANISOEDGE_BASE(I, sigma, iter, max);
%      
%% Acknowledgment
% This function uses the freely available C function developped by
%     <mailto:kranenbu@bigpine.csee.usf.edu Christine Kranenburg>
% available at <http://marathon.csee.usf.edu/edge/>.
%
%% Remark
% The original implementation only thresholded the top and left neighbours. 
% Using 4 neighbors instead of 2 produces thicker but more continuous edges. 
% The non-max then performs edge thinning.
% 
%% See also
% Related:
% <ANISOEDGE.html |ANISOEDGE|>,
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>,
% <CANNYEDGE_BASE.html |CANNYEDGE_BASE|>,
% <CANNYEDGEMAP_BASE.html |CANNYEDGEMAP_BASE|>,
% <ROTHWELLEDGE_BASE.html |ROTHWELLEDGE_BASE|>,
% <CONGRUENCYEDGE_BASE.html |CONGRUENCYEDGE_BASE|>,
% <COMPASSEDGE_BASE.html |COMPASSEDGE_BASE|>,
% <SDGDEDGE_BASE.html |SDGDEDGE_BASE|>,
% <ELDERZUCKEREDGE_BASE.html |ELDERZUCKEREDGE_BASE|>,
% <KOETHEDGE_BASE.html |KOETHEDGE_BASE|>,
% <PETROUEDGE_BASE.html |PETROUEDGE_BASE|>.
% Called:
% ANISOEDGE_MEX.

%% Function implementation
function [edgemap,varargout] = anisoedge_base(I, sigma, iter, max)

%%
% dealing with multispectral images

C = size(I,3);
if C>1
    edgemap = false(size(I));
    if nargout==2,  varargout{1} = zeros(size(I));  end;
    for c=1:C
        [edgemap(:,:,c), tmp] = anisoedge_base(I(:,:,c), sigma, iter, max);
        if nargout==2,  varargout{1}(:,:,c) = tmp;   end;
    end
    return;
end

%%
% call the mex file

[edgemap, mag] = anisoedge_mex(I, sigma, iter, max);
edgemap = ~edgemap;

if nargout==2,  varargout{1} = mag;  end;

end % end of anisoedge_base
