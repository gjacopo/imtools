%% SUSANCORNER - Edge/corner detection or image smoothing.
%
%% Description
% Edge/corner detection or image smoothing using the SUSAN technique of [SB97] 
% already implemented in [SUSAN].
% 
%% Syntax
%     Smap = SUSANCORNER(I);
%     [Smap, Spt] = SUSANCORNER(I, 'c');
%     Smap = SUSANCORNER(I, mode, 'Property', propertyvalue, ...);
% 
%% Inputs
% *|I|* : an input image with size |(X,Y,C)|, where |C>1| when |I| is multichannel.
%
% *|mode|* : optional string for chosing the computing mode, on which the
%     outputs of the SUSAN filter will depend; it is either:
%
% * |'e'| for detecting edges,
% * |'ei'| for edges overlaid on input image,
% * |'s'| for smoothing mode,
% * |'c'| for extracting corners - output a |(n,2)| corners matrix with the 
%          coordinates of the |n| corners in the image,
% * |'ci'| for corners overlaid on the input image;
%
% default: |mode='e'|, ie the edges are extracted.
%
% Property [propertyname  propertyvalues]:
% *|'thres' : optional parameter for brightness threshold, used by all
%     modes; default: |thres=20|;
%
% *|'md'|* : distance mask, either a threshold (with |'s'| mode) or the
%     string |'flat'| for flat |(3,3)| masks (with |'e'| and |'s'| modes).
%
% *|'n'|* : optional boolean for avoiding the post-processing on the binary
%     edge map (runs much faster) with |'e'| modes; default: |n=false|;
%
% *|'q'|* : optional boolean for faster (and usually stabler) |'c'| modes; 
%     edge-like corner suppression not carried out; default: |q=true|.
% 
%% Outputs
% *|Smap|* : depending on the computing mode selected (see above), it can be
%     either:
%
% * a logical (binary) array (mask); this is the default type used with the
%          |'e'| mode,
% * an image of same size as |I| when either the |'s'| mode or the |'ei'|
%          and |'ci'| modes are selected (type |uint8| with overlaid features
%          in the latter cases),
% * a |(n,2)| corners matrix with the coordinates of the |n| corners in the
%          image.
%
% *|Spt|* : coordinates of the corner when |mode='c'|.
%
%% References
% [SB97]  S.M. Smith and J.M. Brady: "{SUSAN} - A new approach to low  
%      level image processing", International Journal of Computer Vision,
%      23(1):45-78, 1997.
%      <http://www.lems.brown.edu/vision/courses/image-processing/Readings/smith95susan.pdf>
%
% [SUSAN]  Source code available at <http://www.fmrib.ox.ac.uk/~steve/susan/>
%
%% See also 
% Related:
% <CORNER.html |CORNER|>,
% <EDGECORNER.html |EDGECORNER|>,
% <HARRISCORNER.html |HARRISCORNER|>,
% <FASTCPDA.html |FASTCPDA|>,
% <FASTCORNER.html |FASTCORNER|>,
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>.
% Called: 
% <SUSANCORNER_BASE.html |SUSANCORNER_BASE|>.

%% Function implementation
function [Smap, varargout] = susancorner(I,varargin)

%% 
% parsing and checking parameters
error(nargchk(1, 18, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

if ~isnumeric(I)
    error('susancorner:inputerror','matrix required in input'); 
end

p = createParser('SUSANCORNER');   
% mandatory parameter
% optional parameters
p.addOptional('mode', 'e', @(x)ischar(x) && ...
    any(strcmpi(x,{'e','ei','c','ci','s'})));
p.addParamValue('thres', 20, @(x)isscalar(x) && x>=0);
p.addParamValue('md', 'flat', @(x)(ischar(x) && strcmp(x,flat)) || ...
    (isscalar(x) && x>0));
p.addParamValue('n', false, @(x)isscalar(x) && (x==true||x==false));
p.addParamValue('q', true, @(x)isscalar(x) && (x==true||x==false));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%% 
% main computation

C = size(I,3);

if strcmpi(p.mode,'c') && nargout==2
    [Smap, varargout{1}] = susancorner_base(I, p.mode, p.thres, p.md, p.n, p.q);
else
     Smap = susancorner_base(I, p.mode, p.thres, p.md, p.n, p.q);
end
    
%%
% display

if p.disp
    figure, 
     if any(strcmpi(p.mode,{'e','¨ei','s'})), 
         imagesc(Smap);
     else
         tmp = false(size(I));
         for c=1:C,  
             tmp(Smap{c}(:,1) + X*(Smap{c}(:,2)-1) + X*Y*(c-1)) = true; 
         end
         imagesc(tmp);
     end
        axis image off;
        if C==1, colormap gray;  end;
        title('corners detected with SUSAN');
end

end % end of susancorner
