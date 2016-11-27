%% MORPHPROFILE - (Derivative) morphological profile by opening and/or closing.
%
%% Description
% Compute the morphological profile, the derivative morphological profile by
% opening and closing (or both) and the morphological multiscale characteristics.
%
%% Syntax
%     DMP = MORPHPROFILE(I);
%     [DMP, Phi] = MORPHPROFILE(I, op, se);
%     [DMP, Phi, se, MP] = MORPHPROFILE(I, op, se, s1[, s2]);
%
%% Input
% *|I|* : input image of size |(X,Y,C)|, with |C>1| when |I| is multispectral.
% 
%% Outputs
% *|MP|* : (derivative of the) morphological profile.
%
% *|Phi|* : morphological multi-scale characteristics.
%
%% References
% [PB00]  M. Pesaresi and J.A. Benediktsson: "Image segmentation based on
%      the derivative of the morphological profile", in Mathematical 
%      Morphology and Its Applications to Image and Signal Processing, 
%      J. Goustsias et al. eds. Norwell (MA, USA), Kluwer, 2000. 
%      <http://www.springerlink.com/content/k217533885154776/>
%
% [PB01]  M. Pesaresi and J.A. Benediktsson: "A new approach for the 
%      morphological segmentation of high-resolution satellite imagery", 
%      IEEE Trans. Geosci. Remote Sens., 39(2):309-320, 2001.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=905239>
%
% [SP02]  P. Soille and M. Pesaresi: "Advances in mathematical morphology 
%      applied to geoscience and remote sensing", IEEE Trans. Geosci. Remote 
%      Sens., 40(9):2042-2055, 2002. 
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1046853>
%
% [Soille03]  P. Soille: "Morphological Image Analysis - Principles and 
%      Applications", 2nd ed. Berlin (Germany), Springer-Verlag, 2003. 
%
% [BPA03]  J.A. Benediktsson, M. Pesaresi and K. Arnason: "Classification
%      and feature extraction for remote sensing images from urban areas based
%      on morphological transformations", IEEE Trans. Geosci. Remote Sens.,
%      41(9):1940-1949, 2003. 
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1232208>
%
%% See also
% Related:
% <matlab:web(whichpath('WATERSHED')) |WATERSHED|>,
% <GRANULOMETRY.html |GRANULOMETRY|>,
% <ASF.html |ASF|>.
% Called:
% <MORPHPROFILE_BASE.html |MORPHPROFILE_BASE|>.

%% Function implementation
function [MP, Phi, varargout] = morphprofile( I, varargin )

%%
% parsing parameters

if isempty(ver('images'))
    error('morphprofile:errortoolbox', 'Image Processing toolbox required');
end

error(nargchk(1, 16, nargin, 'struct'));
error(nargoutchk(1, 4, nargout, 'struct'));


if ~isnumeric(I)
    error('morphprofile:inputparameter','a matrix is required in input'); 
end

p = createParser('MORPHPROFILE');   
p.addOptional('op', 'roc', @(x)ischar(x) && ...
    any(strcmpi(x,{'ro', 'ropen', 'rc', 'rclose', 'roc','rocmax', ...
    'o', 'open', 'c', 'close', 'oc', 'ocmax', ...
    'e', 'erode', 'd', 'dilate', 'ed', 'edmax'})));
p.addOptional('se', 'disk', @(x) (iscell(x) && strcmp(class(x{1}),'strel')) || ...
    (ischar(x) && any(strcmpi(x,{'disk','rectangle','square','diamond', ...
                'line','periodicline','arbitrary','octagon','pair'}))));
p.addParamValue('s1', 2:2:10, @(x)isnumeric(x));
p.addParamValue('s2', [], @(x)isnumeric(x));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% checking/setting variables

C = size(I,3);

if ischar(p.se)
    shape = p.se;
    N = length(p.s1);
    
    if ~iscell(p.s1),
        if any(strcmp(shape,{'disk'}))
            p.s1 = floor(p.s1/2);
        end
        p.s1 = num2cell(p.s1);  
    end
    if ~iscell(p.s2),
        if isempty(p.s2)
            p.s2 = cell(N,1);  
        else
            if isscalar(p.s2),  p.s2 = repmat(p.s2, 1, N);  end;
            p.s2 = num2cell(p.s2); 
        end
    end
    
    p.se = cell(N,1);
    for i=1:N
        p.se{i} = flatstrel(shape, p.s1{i}, p.s2{i});
    end
    
end

%%
% main processing

if nargout==4
    [MP, Phi, varargout{2}] = morphprofile_base( I, p.op, p.se );
else
    [MP, Phi] = morphprofile_base( I, p.op, p.se );
end

if nargout>=3,  varargout{1} = p.se;  end

%%
% display

if p.disp
    
    if any(strcmp(p.op,{'ro','ropen'})),       p.op = 'open-by-rec';
    elseif any(strcmp(p.op,{'rc','rclose'})),  p.op = 'close-by-rec';
    elseif strcmp(p.op,'roc'),                 p.op = 'open- and close-by-rec';
    elseif strcmp(p.op,'rocmax'),              p.op = 'max(open-,close-by-rec)';
    elseif any(strcmp(p.op,{'o','open'})),     p.op = 'open';
    elseif any(strcmp(p.op,{'c','close'})),    p.op = 'close';
    elseif strcmp(p.op,'oc'),                  p.op = 'open-close';
    elseif strcmp(p.op,'ocmax'),               p.op = 'max(open,close)';
    elseif any(strcmp(p.op,{'e','erode'})),    p.op = 'erosion';
    elseif any(strcmp(p.op,{'d''dilate'})),    p.op = 'dilation';
    elseif strcmp(p.op,'ed'),                  p.op = 'erosion-dilation';
    elseif strcmp(p.op,'edmax'),               p.op = 'max(erosion,dilation)';
    end
    
    N = length(MP);
    if any(strcmp(p.op,{'roc','oc','ed'})),  N = N/2;  end
    figure; ndisp = 3; mdisp = ceil(N/ndisp);
    for i=1:N
        subplot(mdisp, ndisp, i),  imagesc(rescale(MP{i})), axis image off;
    end
    if C==1,  colormap gray; end
    if any(strcmp(p.op,{'roc','oc','ed'}))
        figure
        for i=1:N
            subplot(mdisp, ndisp, i),  imagesc(rescale(MP{2*N-i+1})), axis image off;
        end
        if C==1,  colormap gray; end
    else
        suptitle(['DMP of ' p.op]);
    end
    figure, imagesc(rescale(Phi)), axis image off, colormap gray,
    title('morphological multi-scale characteristics');
end

end % end of morphprofile 
