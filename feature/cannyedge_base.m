%% CANNYEDGE_BASE - Base function for CANNYEDGE.
% 
%% Syntax:
%     edgemap = CANNYEDGE_BASE(I, sigma, der, samp, hyst, reduce);
%     [edgemap, mag, orient] = CANNYEDGE_BASE(I, sigma, der, samp, hyst, reduce);
%
%% See also
% Related:
% <CANNYEDGE.html |CANNYEDGE|>,
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>,
% <ROTHWELLEDGE_BASE.html |ROTHWELLEDGE_BASE|>,
% <CONGRUENCYEDGE_BASE.html |CONGRUENCYEDGE_BASE|>,
% <COMPASSEDGE_BASE.html |COMPASSEDGE_BASE|>,
% <ANISOEDGE_BASE.html |ANISOEDGE_BASE|>,
% <ELDERZUCKEREDGE_BASE.html |ELDERZUCKEREDGE_BASE|>,
% <KOETHEDGE_BASE.html |KOETHEDGE_BASE|>,
% <SDGDEDGE_BASE.html |SDGDEDGE_BASE|>,
% <PETROUEDGE_BASE.html |PETROUEDGE_BASE|>.
% Called:
% <matlab:webpub(whichpath('EDGE')) |EDGE|>,
% <matlab:webpub(whichpath('BWMORPH')) |BWMORPH|>,
% <matlab:webpub(whichpath('RGB2GRAY')) |RGB2GRAY|>,
% <../../vista/html/CANNYEDGES.html |CANNYEDGES|>,
% <../../derive/html/GRDSMOOTH_BASE.html |GRDSMOOTH_BASE|>,
% <../../geometry/html/UPSCALEXY_BASE.html |UPSCALEXY_BASE|>. 

%% Function implementation
function [edgemap,varargout] = cannyedge_base(I, sigma, der, samp, hyst, reduce)

%% 
% checking/setting variables
[X,Y,C] = size(I);

%% 
% dealing with multispectral images

if C>1 
    if (islogical(reduce) && ~reduce) || ...
            (ischar(reduce) && strcmpi(reduce,'eor'))
        % process channel by channel
        edgemap = false(X*samp,Y*samp,C);
        for i=1:nargout-1
            varargout{i} = zeros(X*samp,Y*samp,C);
        end
        for c=1:C
            [edgemap(:,:,c), tmp1, tmp2] = ...
                cannyedge_base(I(:,:,c), sigma, der, samp, hyst, reduce);
            if nargout>=2, varargout{1}(:,:,c) = tmp1; end
            if nargout==3, varargout{2}(:,:,c) = tmp2; end
        end
        if strcmpi(reduce,'eor')
            for c=1:C
                edgemap(:,:,1) = edgemap(:,:,1) | edgemap(:,:,c);
            end
            edgemap = edgemap(:,:,1);
            if ~isempty(ver('images'))
                edgemap = bwmorph(edgemap,'thin',Inf);
            end
        end
        return;
        
    elseif islogical(reduce) || ...
            (ischar(reduce) && any(strcmpi(reduce,{'isum','igray','imax'})))
        % at that point, if reduce is logical, it is necessarly true
        
        if C==3 && strcmpi(reduce,'igray')
            if ~isempty(ver('images')),  I = rgb2gray(I);
            else    I = 0.2989*I(:,:,1) + 0.587*I(:,:,2) + 0.114*I(:,:,3);
            end
            I = rescale(I,0,1); 
           
        elseif islogical(reduce) || strcmpi(reduce,'isum')
            I = rescale(sum(I,3),0,1);
            
        elseif strcmpi(reduce,'imax')
            I = rescale(max(I,[],3),0,1);
        end
        
    end
    
end

%% 
% computing edges

if samp>1
    S =[samp samp];
    % perform interpolation channel by channel
    I = upscalexy_base(I, S, 'linear');
end

if any(strcmpi(der,{'matlab','vista','edge'})) && ...
        ~(ischar(reduce) && strcmpi(reduce,'gmax'))
    
    % already implemented edge detection
   if any(strcmpi(der,{'matlab','edge'}))
        edgemap = edge(I, 'canny', [], sigma);      
        mag = []; or = [];
    
   elseif strcmpi(der,'vista')   
       [edgemap, or, mag] = cannyedges(I, sigma);
   end
   
else % all other cases
    
    % if strcmpi(der,'kovesi'),  der = 'fleck';  end;
    [gx, gy] = grdsmooth_base(I, sigma, der, [], 'xy');

    if C>1 && ischar(reduce) && strcmpi(reduce,'gmax');
        gx = max(gx, [], 3);
        gy = max(gy, [], 3);
    end
    
    if strcmpi(der,'fleck'),  dermap = 'kovesi';
    else                      dermap = der;
    end
    [edgemap, mag, or] = cannyedgemap_base(gx, gy, dermap, [], [], hyst, []);
    
end

%%
% special case when multispectral

if C>1 && ischar(reduce) && strcmpi(reduce,'eor')
    [edgemap, J] = max(edgemap, [], 3);
    subm = reshape(1:X*Y,[X Y]) + (J-1)*X*Y;
    mag = mag(subm);
    or = or(subm);
end

%%
% specify outputs

if nargout>=2,     varargout{1} = mag;     end
if nargout==3,     varargout{2} = or;     end

end % end of cannyedge_base
