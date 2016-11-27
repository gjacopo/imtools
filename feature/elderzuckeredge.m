%% ELDERZUCKEREDGE - Multi-scale edge detector.
%
%% Description
% This function runs the Elder and Zucker multi-scale edge detector of [EZ98] 
% and returns the edgel image.
%
%% Syntax
%     edgemap = ELDERZUCKEREDGE(I);  
%     [edgemap,scmap,blur] = ELDERZUCKEREDGE(I, sigma);  
%
%% Inputs
% *|I|* : input image.
%      
% *|sigma|* : an estimate of the sensor noise of the image.
%
%% Outputs
% *|edgemap|* :  an edge map, not necessarily optimized.
%      
% *|scmap|* : optional output storing the minimum reliable scale of the second
%     derivative estimator.
%      
% *|blur|* : optional output storing the blur estimate at every pixel.
%
%% References
% [EZ98]  J.H. Elder and S.W. Zucker: "Local scale control for edge 
%      detection and blur estimation", IEEE Trans. on Pattern Analysis and
%      Machine Intelligence, 20(7), 1998.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=689301&tag=1>
%      
% [Elder99]  J.H. Elder: "Are edges incomplete?", International Journal of
%      Computer Vision, 1999. 
%      <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.27.9846>
%
%% See also
% Related: 
% <EDGECORNER.html |EDGECORNER|>,
% <CANNYEDGE.html |CANNYEDGE|>,
% <CANNYEDGEPROD.html |CANNYEDGEPROD|>,
% <CANNYEDGEMAP.html |CANNYEDGEMAP|>,
% <SDGDEDGE.html |SDGDEDGE|>,
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>,
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <ANISOEDGE.html |ANISOEDGE|>,
% <ROTHWELLEDGE.html |ROTHWELLEDGE|>,
% <KOETHEDGE.html |KOETHEDGE|>,
% <PETROUEDGE.html |PETROUEDGE|>.
% Called: 
% <ELDERZUCKEREDGE_BASE.html |ELDERZUCKEREDGE_BASE|>.

%% Function implementation
function [edgemap, varargout] = elderzuckeredge(I, varargin)

%%
% parsing parameters

error(nargchk(1, 12, nargin, 'struct'));
error(nargoutchk(1, 3, nargout, 'struct'));

if ~isnumeric(I)
    error('elderzuckeredge:inputparameter','a matrix is required in input'); 
end

p = createParser('ELDERZUCKEREDGE');   
p.addOptional('sigma', 0.05, @(x)isscalar(x) && x>=0.05);
p.addParamValue('reduce', false, @(x)islogical(x) || ...
    (ischar(x) && any(strcmpi(x,{'igray','imax','isum','gmax','eor'}))));
% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% main calculation
if nargout==1
    edgemap = elderzuckeredge_base(I, p.sigma, p.reduce);

else
    [edgemap, scmap, blur] = elderzuckeredge_base(I, p.sigma, p.reduce);    
    if nargout>=2,     varargout{1} = scmap;     end
    if nargout==3,     varargout{2} = blur;     end
end

%%
% display

if p.disp
    figure, imagesc(edgemap), axis image off, title('EZ edge map');
    if size(edgemap,3) == 1, colormap 'gray'; end;
    if nargout==3
        figure, subplot(1,2,1)
        imshow(scmap), axis image off, title('minimum reliable scale');
        subplot(1,2,2)
        imshow(rescale(blur)), axis image off, title('blur estimate');
    end
end

end % end of elderzuckeredge