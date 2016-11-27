%% WMEDFILT2 - Weighted median filtering.
% 
%% Description
% Perform weighted median filtering of a (possibly multispectral) image.
%
%% Syntax
%     WMI = WMEDFILT2(I,wk);
%
%% Inputs
% *|I|* : input image, possibly multispectral.
%
% *|wk|* : kernel defining the neighbourhood and the local weights used for 
%      weighted median estimation; a negative weight indicates the position
%      of the central pixel tested, otherwise the kernel is supposed to be
%      square-shaped with odd size and the barycenter is considered.
%      
%% Output
% *|WMI|* : weighted median filtered image with same size as the input |I|
%      each pixel in |WMI| contains the median value in the neighbourhood of
%      |I| defined by the non-negative values of |wk| around the corresponding
%      pixel in |I|. 
%
%% References
% [Brown84]  D.R.K. Brownrigg: "The weighted median filter", Image 
%      Processing and Computer Vision, Communications of the ACM, 27(8):
%      807-818, 1984.
%      <http://dl.acm.org/citation.cfm?id=358222&bnc=1>
%
% [Yaros94]  L. Yaroslavsky: "Local criteria: a unified approach to local 
%      adaptive linear and rank filters for image restoration and enhancement", 
%      Proc. IEEE ICIP, pp. 517-521, 1994.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=413624&tag=1>
%
% [Yaros96]  L. Yaroslavsky: "Local adaptative filters for image 
%      restauration and enhancement", LNCIS, vol. 219, pp. 31-39, 1996.
%      <http://www.springerlink.com/content/y45n7228700540n2/>
% 
%% See also
% Related:
% <matlab:webpub(whichpath('MEDIAN')) |MEDIAN|>,
% <matlab:webpub(whichpath('MEDFILT2')) |MEDFILT2|>,
% <matlab:webpub(whichpath('ORDFILT2')) |ORDFILT2|>,
% <WORDFILT2.html |WORDFILT2|>.
% Called:
% <WMEDFILT2_BASE.html |WMEDFILT2_BASE|>.

%% Function implementation
function WMI = wmedfilt2(I,wk)

error(nargchk(2, 2, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

%% 
% parsing and checking parameters
if ~isnumeric(I)
    error('wmedfilt2:inputerror','matrice is required in input');
end

p = createParser('WMEDFILT2');   % create an instance of the inputParser class.
p.addRequired('wk', @(x) isequal(round(x),x));

% parse and validate all input arguments
p.parse(wk); 
p = getvarParser(p);                                                            

%%
% internal variables and further checking

C = size(I,3);
[x y c] = size(p.wk);

if ~isequal(round(p.wk),p.wk)
    warning('wmedfilt2:inputwarning',...
        'the input weighting matrix must have integer values');
end

if c~=1 && c~=C
    error('wmedfilt2:inputargument', ...
    'input kernel must have 1 channel or the same number as the input');
end

cpos = find(p.wk<0,1,'first');
if isempty(cpos) && (x~=y || mod(x,2)~=1) 
    error('wmedfilt2:inputargument',...
        'input kernel must be square shaped with odd size');
end

%% 
% main processing

WMI = wmedfilt2_base(I, p.wk);

%%
% display

if p.disp
    figure, imagesc(rescale(WMI,0,1)), axis image off;
    if C==1,  colormap gray; end;
    title('weighted median filtered image');
end

end % end of wmedfilt2
