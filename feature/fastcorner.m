%% FASTCORNER - Wrapper of functions for fast corner detection.
%
%% Description
% Nothing else than a wrapper for the functions [FAST] of the FAST method of
% [RD05,RD06] used for corner detection.
%
%% Syntax
%    ptcorner = FASTCORNER(I);
%    [ptcorner, cornermap] = FASTCORNER(I, 'Property', propertyvalue,... );
%
%% Inputs
% *|I|* : image.
%      
% *|numa|* : (optional) string or scalar value setting the fast method to be
%     used; it is any number/char in |{9,10,11,12}|; default: |numa='9'|.
%      
% *|nonmax|* : (optional) logical flag for operating non maximum suppression;
%     default: |nonmax=true|.
%      
% *|thres|* : thres to apply; default: |thres=1|.
%
%% Outputs
% *|ptcorner|* : list of corner points.
%      
% *|cornermap|* : 2D spatial representation of ptcorner; cornermap is a logical
%     matrix.
%
%% References
% [RD05]  E. Rosten and T. Drummond: "Fusing points and lines for high 
%      performance tracking", Proc. ICCV, vol. 2, pp 1508-1511, 2005.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1544896>
%      
% [RD06]  E. Rosten and T. Drummond: "Machine learning for high-speed
%      corner detection", Proc. ECCV, vol. 1, pp 430-443, 2006.
%      
% [FAST]  Source code available at <http://mi.eng.cam.ac.uk/~er258/work/fast.html>
%      and <http://www.mathworks.com/matlabcentral/fileexchange/13006-fast-corner-detector>.
%
%% See also
% Related:
% <CORNER.html |CORNER|>,
% <EDGECORNER.html |EDGECORNER|>,
% <SUSANCORNER.html |SUSANCORNER|>,
% <HARRISCORNER.html |HARRISCORNER|>,
% <FASTCPDA.html |FASTCPDA|>,
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>.
% Called: 
% <FASTCORNER_BASE.html |FASTCORNER_BASE|>.
    
%% Function implementation
function [ptcorner,cornermap] = fastcorner(I,varargin)

% %% Check if possible
% if ~exist('fast_mex','file')
%     error('mex file fast_mex not found');
% end

%%
% parsing parameters

error(nargchk(1, 12, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(I)
    error('fastcorner:inputerror','matrix required in input'); 
end

% optional parameters
p = createParser('FASTCORNER');   
% principal optional parameters
p.addOptional('numa', 9, @(x) (isscalar(x) && ismember(x,[9,10,11,12])) || ...
    (ischar(x) && any(strcmpi(x,{'9','10','11','12'}))));
p.addOptional('nonmax', true, @(x)islogical(x));
% additional optional parameters
p.addOptional('thres', 1, @(x)isscalar(x) && x>=0);

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%% 
% setting internal variables

if ischar(p.numa), p.numa = str2double(p.numa);  end;
    
%%
% main processing

[ptcorner, cornermap] = fastcorner_base(I, p.numa, p.nonmax, p.thres);

%%
% display

if p.disp
    figure, imagesc(cornermap), axis image off;
    if size(cornermap,3)==1, colormap gray;  end;
    title('corners detected with FAST');
end

end % end of fastcorner