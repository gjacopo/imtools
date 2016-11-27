%% ROTHWELLEDGE_BASE - Base function for ROTHWELLEDGE.
%
%% Syntax
%    [edgemap, mag, or] = ROTHWELLEDGE_BASE(I, sigma, low, alpha, samp, reduce);
%    [edgemap, mag, or] = ROTHWELLEDGE_BASE(gx, gy, low, alpha, samp, reduce);
%      
%% Acknowledgment
% This function uses the C function developped Heath et al. for the comparative
% study in [HSSB97]; the mex file in |ROTHWELL_MEX| calls directly the code
% made available by the authors in the page: 
% <ftp://figment.csee.usf.edu/pub/Edge_Comparison/source_code/rothwell.src>
%
%% See also
% Related:
% <ROTHWELLEDGE.html |ROTHWELLEDGE|>,
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>,
% <CANNYEDGE_BASE.html |CANNYEDGE_BASE|>,
% <CONGRUENCYEDGE_BASE.html |CONGRUENCYEDGE_BASE|>,
% <COMPASSEDGE_BASE.html |COMPASSEDGE_BASE|>,
% <ANISOEDGE_BASE.html |ANISOEDGE_BASE|>,
% <ELDERZUCKEREDGE_BASE.html |ELDERZUCKEREDGE_BASE|>,
% <KOETHEDGE_BASE.html |KOETHEDGE_BASE|>,
% <SDGDEDGE_BASE.html |SDGDEDGE_BASE|>,
% <PETROUEDGE_BASE.html |PETROUEDGE_BASE|>.
% Called: 
% ROTHWELLEDGE_MEX.
 
%% Function implementation
function [edgemap,varargout] = rothwelledge_base(I, V, low, alpha, samp, reduce)

if nargin<5,  samp = 1;  
    if nargin<6,   reduce = false;
    end
end; % we allow variable number of inputs

%% 
% dealing with multispectral images

[X,Y,C] = size(I);

if C>1
    if (islogical(reduce) && ~reduce) || ...
            (ischar(reduce) && strcmpi(reduce,'eor'))
        % process channel by channel
        edgemap = false(X*samp, Y*samp, C);
        for i=1:nargout-1,
            varargout{i} = zeros(X*samp, Y*samp, C);
        end
        if isscalar(V), V = repmat(V, [1 1 3]);  end
        for c=1:C
            [edgemap(:,:,c), mag, or] = ...
                rothwelledge_base(I(:,:,c), V(:,:,c), low, alpha, samp);
            if nargout>=2,  varargout{1}(:,:,c) = or;   end
            if nargout==3,  varargout{2}(:,:,c) = mag;  end
        end
        if strcmpi(reduce,'eor')
            for c=1:C
                edgemap(:,:,1) = edgemap(:,:,1) | edgemap(:,:,c);
            end
            edgemap = bwmorph(edgemap(:,:,1),'thin',Inf);
        end
        return;
        
    elseif islogical(reduce) || ...
            (ischar(reduce) && any(strcmpi(reduce,{'isum','igray','imax'})))
        % at that point, if reduce is logical, it is necessarly true
        
        if C==3 && strcmpi(reduce,'igray')
            I = rgb2gray(rescale(I,0,1));
            
        elseif islogical(reduce) || strcmpi(reduce,'isum')
            if C==3
                I = 0.29*I(:,:,1) + 0.59*I(:,:,2) + 0.11*I(:,:,3);
                I = rescale(I,0,1);
            else
                I = sum(rescale(I,0,1), 3);
            end
            
        elseif strcmpi(reduce,'imax')
            I = max(rescale(I,0,1), [], 3);
        end
        
    end
    
end

I = rescale(I,0,255); % requirement when calling rothwelledge_mex...

%% 
% set some internal parameter fo padding the input image (see function

if samp>1
    S =[samp samp];
    % perform interpolation channel by channel
    I = upscalexy_base(I, S, 'linear');
end

%%
% |smooth_image| of |rothwell_mex|
if isscalar(V)
    sigma = samp * V;
    gy = [];
    gauss_tail = 0.015;
    kwidth = ceil(sigma*sqrt(2*log(1/gauss_tail))+1); %
    A = padarray(I, [kwidth kwidth], 'symmetric', 'both');
    
elseif isnumeric(V)
    sigma = [];
    gy = V;
    if samp>1,  gy = upscalexy_base(gy, S, 'linear');  end
    A = I;
    kwidth = 0;
end

%% 
% call the mex file

[edgemap, or, mag] = rothwelledge_mex(A, gy, sigma, low, alpha);

%%
% crop
if isscalar(V)
    edgemap = edgemap(kwidth+1:kwidth+X*samp,kwidth+1:kwidth+Y*samp);
end

%%
% outputs
if nargout>=2
    varargout{1} = mag(kwidth+1:kwidth+X*samp,kwidth+1:kwidth+Y*samp);
end
if nargout==3 
    varargout{2} = or(kwidth+1:kwidth+X*samp,kwidth+1:kwidth+Y*samp);
end

end % end of rothwelledge_base