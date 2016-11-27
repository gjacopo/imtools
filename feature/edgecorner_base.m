%% EDGECORNER_BASE - Base function for EDGECORNER.
% 
%% Syntax
%     [edgemap, cormap] = ...
%              EDGECORNER_BASE(I, edge, corner, rho, sigma, thres, reduce);
%     [edgemap, cormap, MO] = ...
%              EDGECORNER_BASE(I, edge, corner, rho, sigma, thres, reduce);
%     [edgemap, cormap, T] = ...
%              EDGECORNER_BASE(I, edge, corner, rho, sigma, thres, reduce);
%
%% See also
% Related: 
% <EDGECORNER.html |EDGECORNER|>,
% <CORNER_BASE.html |CORNER_BASE|>.
% Called: 
% <matlab:webpub(whichpath('EDGE')) |EDGE|>,
% <../../derive/html/GRDSMOOTH_BASE.html |GRDSMOOTH_BASE|>,
% <../../algebra/html/RESCALE.html |RESCALE|>. 
% <../../../vista/html/CANNYEDGES.html |CANNYEDGES|>,
% <CANNYEDGE_BASE.html |CANNYEDGE_BASE|>,
% <ROTHWELLEDGE_BASE.html |ROTHWELLEDGE_BASE|>,
% <CANNYEDGEMAP_BASE.html |CANNYEDGEMAP_BASE|>,
% <CONGRUENCYEDGE_BASE.html |CONGRUENCYEDGE_BASE|>,
% <COMPASSEDGE_BASE.html |COMPASSEDGE_BASE|>,
% <ANISOEDGE_BASE.html |ANISOEDGE_BASE|>,
% <ELDERZUCKEREDGE_BASE.html |ELDERZUCKEREDGE_BASE|>,
% <KOETHEDGE_BASE.html |KOETHEDGE_BASE|>,
% <SDGDEDGE_BASE.html |SDGDEDGE_BASE|>,
% <PETROUEDGE_BASE.html |PETROUEDGE_BASE|>,
% <HARRISCORNER.html |HARRISCORNER|>,
% <SUSANCORNER_BASE.html |SUSANCORNER_BASE|>,
% <FASTCPDA_BASE.html |FASTCPDA_BASE|>,
% <FASTCORNER_BASE.html |FASTCORNER_BASE|>,
% <matlab:webpub(whichpath('BWMORPH')) |BWMORPH|>.

%% Function implementation
function [edgemap, cormap, varargout]  = ...
    edgecorner_base(I, met_edge, met_corner, rho, sigma, thres, red, varargin)

%% 
% checking/setting variables

error(nargchk(1, 9, nargin, 'struct'));
error(nargoutchk(1, 3, nargout, 'struct'));

% note: the varargin are used with the method met_edge='koethe' only
if nargin>=8,  int = varargin{1};  else   int = 'ani'; end
if nargin==9,  samp = varargin{2};  else samp=2;     end
    
%%
% we foresee that edgemap or cormap may be empty in the case the 'met_'
% variables are set to false

%%
% dealing with multispectral images

C = size(I,3);

if C>1 && ...
        ( (~strcmpi(met_edge,'koethe') && ...
        ((islogical(red)&&~red) || (ischar(red)&&strcmpi(red,'eor')))))
    if ischar(met_edge),  edgemap = zeros(size(I)); % else set to dummy
    else                  edgemap = zeros(1,1,3);  end
    if ischar(met_corner),  cormap = zeros(size(I)); % else set to dummy
    else                    cormap = zeros(1,1,3);  end
    if nargout==3,  varargout{1} = cell(C,1);  end
    for c=1:C
        [A, B, tmp] = ...
            edgecorner_base(I(:,:,c), met_edge, met_corner, rho, sigma, ...
            thres, red, int, samp);
        if ischar(met_edge),  edgemap(:,:,c) = A;  end
        if ischar(met_corner),  cormap(:,:,c) = B; end
        if nargout==3,  varargout{1}{c} = tmp; end
    end
    if strcmpi(red,'eor')
        for c=2:C
            edgemap(:,:,1) = edgemap(:,:,1) | edgemap(:,:,c);
            cormap(:,:,1) = cormap(:,:,1) | cormap(:,:,c);
            varargout{1}{1}(:,:,1) = ...
                max(cat(3,varargout{1}{1}(:,:,1),varargout{1}{c}(:,:,1)),[],3);
            varargout{1}{1}(:,:,2) = ...
                varargout{1}{1}(:,:,2) + varargout{1}{c}(:,:,2);
        end
        varargout{1}{1}(:,:,2) = varargout{1}{1}(:,:,2) / 3;
        varargout{1} = varargout{1}{1}; % reduce the cell to a single matrix
        edgemap = edgemap(:,:,1);
        if ~isempty(ver('images'))  
            edgemap = bwmorph(edgemap,'thin',Inf);
        end
        cormap = cormap(:,:,1);
    end
    if ~ischar(met_edge),  edgemap = []; end
    if ~ischar(met_corner),  cormap = []; end
    return;
end

if C>1 && ischar(red) %&& any(strcmpi(red,{'isum','igray','imax'}))
    
    if C==3 && strcmpi(red,'igray')
            I = rgb2gray(rescale(I,0,1));
    
    elseif strcmpi(red,'isum') 
        if C==3
            I = 0.3*I(:,:,1) + 0.59*I(:,:,2) + 0.11*I(:,:,3);
            I = rescale(I,0,1);
        else
           I = sum(rescale(I,0,1), 3);
        end
       
    elseif strcmpi(red,'imax')
        I = max(rescale(I,0,1), [], 3);
    end

end

%%
% computing edges and junctions all at once: congruency, compass and Koethe
% methods

if strcmpi(met_edge,'congrue')
    [M, m, or]  = ...
        congruencyedge_base(I, 4, 6, 3, 2.1, 0.55, 2, 0.5, 10, -1);
    if numel(thres)==1,  thres = [thres, thres];  end
    edgemap = M>thres(1);
    cormap = m>thres(2); 
    mag = [];
    
elseif strcmpi(met_edge,'compass')   
    [S, O]  = compassedge_base(I, sigma, 180, 6, (size(I,3)~=3 || red));  %#ok
    edgemap = S>thres(1);  
    % use PCCORNER
    mag = []; or = [];
    
elseif strcmpi(met_edge,'koethe')
    % note : samp=2 and int='ani' corresponds to the method proposed by
    % Koethe
    [edgemap, cormap, mag] = ...
        koethedge_base(I, rho, sigma, 'fast', int, samp, 'koe', [], 3 );
    % note that in fact, mag is a tensor, it is nothing else than the GST
    or = [];
end

%% 
% computing edges

if any(strcmpi(met_edge,{'canny','log','rothwell','black','elder'})) ||...
        (strcmpi(met_edge,'vista') && ~strcmpi(red,'gmax'))
    
    % already implemented edge detection
    if any(strcmpi(met_edge,{'canny','log'}))
        edgemap = edge(I,met_edge, [], sigma);
        mag = []; or = [];
        
    elseif strcmpi(met_edge,'vista')
        [edgemap, or, mag] = cannyedges(I, sigma);
        
    elseif strcmpi(met_edge,'rothwell')
        [edgemap, mag, or] = rothwelledge_base(I, sigma, 5, 0.8);
        
    elseif strcmpi(met_edge,'black')
        [edgemap, mag] = anisoedge_base(I, sigma, 100, true);
        or = [];
        
    elseif strcmpi(met_edge,'elder')
        edgemap = elderzuckeredge_base(I, sigma);
        mag = [];  or = [];
        
    end
    
    
elseif ischar(met_edge) && ...
        ~any(strcmpi(met_edge,{'congrue','compass','koethe'})) % all other cases
    
    [gx, gy] = grdsmooth_base(I, sigma, met_edge, [], 'xy');
    
    if ischar(red) && strcmpi(red,'gmax');
        gx = max(gx, [], 3);
        gy = max(gy, [], 3);
    end
    
    [edgemap, mag, or] = ...
        cannyedgemap_base(gx, gy, 'matlab', [], [], [], [1/3 0.08]);
    
elseif islogical(met_edge) && ~met_edge
        edgemap = [];
        mag = []; or = [];
end

%% 
% compute junctions if not already done

if ischar(met_corner) && ...
        (~exist('cormap','var') || isempty(cormap))
    if any(strcmp(met_corner,{'harris','noble'})) && exist('gx','var')
        cormap = ...
            corner_base(gx, met_corner, [], 0.06, 3, [], [], [], gy, rho, false);
    else
        cormap = ...
            corner_base(I, met_corner, thres, 0.06, 3, true, 1, 157, sigma, rho, false);
    end
    
elseif islogical(met_corner) && ~met_corner
    cormap = [];

end

%% 
% other possible outputs: we need to compute mag and/or or if not already done

if nargout==3
    if (C==1 || ~strcmpi(met_edge,'koethe')) && (isempty(mag) || isempty(or))
        if exist('gy','var') && exist('gy','var')
            if isempty(mag)
                mag = hypot(gx,gy);
                mag = mag / max(mag(:));
            end
            if isempty(or)
                or = atan2(gy, gx);        
               % or = or.*(or>=0) + (or+pi).*(or<0); % in [0,pi]
            end
        else
            [~, ~, mag1, or1] = grdsmooth_base(I, sigma, 'fast', [], 'xy');
            if isempty(mag),  mag = mag1;  end
            if isempty(or),  or = or1;  end
        end
    end
    if (C==1 || ~strcmpi(met_edge,'koethe')),  mag = cat(3, mag, or);  end
    varargout{1} = mag;
end

end % end of edgecorner_base
