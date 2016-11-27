%% FASTCORNER_BASE - Base function for FASTCORNER.
%
%% Syntax
%     [ptcorner, cornermap] = FASTCORNER_BASE(I, numa, flnonmax, thres);
%
%% Acknowledgment
% This function calls the original FAST functions of E.Rosten distributed at
% <http://mi.eng.cam.ac.uk/~er258/work/fast.html>.
%
%% See also
% Related:
% <FASTCORNER.html |FASTCORNER|>,
% <FASTCPDA_BASE.html |FASTCPDA_BASE|>,
% <HARRISCORNER_BASE.html |HARRISCORNER_BASE|>,
% <SUSANCORNER_BASE.html |SUSANCORNER_BASE|>,
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>.
% Called: 
% FAST_MEX,
% <../../../toolbox/feature/fast/html/FAST_NONMAX.html |FAST_NONMAX|>,
% <../../../toolbox/feature/fast/html/FAST9.html |FAST9|>,
% <../../../toolbox/feature/fast/html/FAST10.html |FAST10|>,
% <../../../toolbox/feature/fast/html/FAST11.html |FAST11|>,
% <../../../toolbox/feature/fast/html/FAST12.html |FAST12|>.

%% Function implementation
function [ptcorner, varargout] = fastcorner_base(I, numa, flnonmax, thres)

%% 
% dealing with multispectral images
[X,Y,C] = size(I);

ptcorner = cell(C,1);

if C>1
    if nargout==2,  varargout{1} = false(size(I));  end;
    for c=1:C
        [tmp1,tmp2] = fastcorner_base(I(:,:,c), numa, flnonmax, thres);
        ptcorner{c} = tmp1{1};
        if nargout==2,  varargout{1}(:,:,c) = tmp2;  end;
    end
    return;
end

%% 
% call the mex file or the matlab function

if exist('fast_mex','file')  && strcmpi(class(I),'uint8') 
    % preferably the mex file
    if ischar(numa),  numa = str2double(numa);  end
    [ptcorner{1}, cornermap] = fast_mex(I, numa, flnonmax, thres);
    
else
    ffast = str2func(['fast' num2str(numa)]);
    ptcorner{1} = ffast(I, thres, flnonmax);
    if nargout==2
        cornermap = false(X,Y);
        cornermap(sub2ind([X,Y],ptcorner{1}(:,2),ptcorner{1}(:,1))) = true;
    end
end

if nargout==2, varargout{1} = cornermap;  end;

end % end of fastcorner_base