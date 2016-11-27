%% SUSANCORNER_BASE - Base function for SUSANCORNER.
% 
%% Syntax
%     Smap = SUSANCORNER_BASE(I, mode, thres, md, n, q);
%     [ptcorner, cornermap] = SUSANCORNER_BASE(I, mode, thres, md, n, q);
%
%% Acknowledgment
% The original C code and algorithm description at:
%    <http://www.fmrib.ox.ac.uk/~steve/susan/>
% SUSAN Version 2l (C) 1995-1997: <mailto:steve@fmrib.ox.ac.uk Stephen Smith>, DRA UK.
% author of the mexification (see susan_mex): , <mailto:jluis@ualg.pt Joaquim Luis>
%
%% See also
% Related:
% <SUSANCORNER.html |SUSANCORNER|>,
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>,
% <HARRISCORNER_BASE.html |HARRISCORNER_BASE|>,
% <FASTCPDA_BASE.html |FASTCPDA_BASE|>,
% <FASTCORNER_BASE.html |FASTCORNER_BASE|>.
% Called: 
% SUSAN_MEX.

%% Function implementation
function [Smap, varargout] = susancorner_base(I, mode, thres, md, n, q)

%%
% dealing with multispectral images

[X,Y,C] = size(I);

if any(strcmpi(mode,{'e','¨ei','s'})),   Smap = false(size(I));
else
    Smap = cell(C,1);  
    if nargout==2,  varargout{1} = false(X,Y);  end;
end

if C>1   
    for c=1:C
        [tmp1, tmp2] = susancorner_base(I(:,:,c), mode, thres, md, n, q);
        if any(strcmpi(mode,{'e','¨ei','s'})),     Smap(:,:,c) = tmp1;
        else
            Smap{c} = tmp1{1};
            if nargout==2,  varargout{1}(:,:,c) = tmp2;  end;
        end;
    end
    return;
end

%% 
% main computation

% susan_mex is available
if exist('susan_mex','file') && strcmpi(class(I),'uint8') %&& false

    % create the string storing the list of arguments
    larg = [ ',''-' num2str(mode) ''',''-t' num2str(thres) ''''];
    if ischar(md) % && strcmp(mask,'flat')
        larg = [larg ',''-3'''];
    elseif any(strcmp(mode,{'e','ei','s'}))
        larg = [larg ',''-d' num2str(md) ''''];
    end
    if n,        larg = [larg ',''-n'''];    end
    if q,        larg = [larg ',''-q'''];    end

    % run the mex file
    eval(['S=susan_mex(I' larg ');']);
    if any(strcmp(mode,{'e','ei','s'}))
        Smap = S;
    else
        Smap{1} = S;
    end
    
elseif any(strcmp(mode,{'e','ei'}))
    % otherwise for matlab: an attempt to do withtout the C
    
    % mask for selecting the pixels within the circular region (37 pixels), as
    % used in the SUSAN algorithm
    mask = [ 0 0 1 1 1 0 0 ,...
        0 1 1 1 1 1 0,...
        1 1 1 1 1 1 1,...
        1 1 1 1 1 1 1,...
        1 1 1 1 1 1 1,...
        0 1 1 1 1 1 0,...
        0 0 1 1 1 0 0];
    mask = mask(:);
    wmask = 3;
    
    % define the USAN area
    nmax = 3*37/4;
    
    % padding the image
    pad = 2*wmask+1;
    A = padarray(I, [wmask wmask],'replicate','both');
    pixA = reshape(1:numel(A),size(A));
    pixI = reshape(pixA(wmask+1:wmask+X,wmask+1:wmask+Y),1,X*Y);
    
    % index of the centered neighbour window of analysis
    indI = -wmask:wmask;
    for i=1:wmask
        indI = [ (-i*X-wmask):(-i*X+wmask), ...
            indI, ...
            (i*X-wmask):(i*X+wmask) ];                                 %#ok
    end
    indI = indI'; indI = indI(:);
    
    indc = wmask * (pad+1) +1;
    
    % the output image indicating found edges
    Smap = zeros(size(A));
    
    for in=pixI
        c = mask .* A(indI + in);
        
        % thresholding scheme: c = exp(-{(I(r)-I(r0))/t}^(5/6)} applied to
        % the current neighborhood of the center pixel within the circle defined
        % by the mask
        tmp = (c-c(indc))/thres;
        tmp = tmp.^6;
        c = exp(-tmp);
        % if binary thresholding is applied
        % c(abs(c-c(indc))>threshold)=0;
        % c(abs(c-c(indc))<=threshold)=1;
        
        g = sum(c(:));
        
        if nmax<g, Smap(in) = g-nmax; end
    end
    
    Smap = Smap(pixI);
    
else
    error('susancorner_base:methoderror', ...
        ['method ' mode ' implemented only using susan_mex']);
    
end

if nargout==2
    if ~strcmpi(mode,'c'),     varargout{1} = [];
    else
        varargout{1} = false(X,Y);
        varargout{1}(sub2ind([X,Y],Smap{1}(:,1),Smap{1}(:,2))) = true;
    end;
end

end % end of susancorner_base
