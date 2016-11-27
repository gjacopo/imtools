%% WORDFILT2 - Weighted order-statistic filtering.
% 
%% Description
% Perform (integer) weighted order-statistic filtering of an image (possibly
% multispectral).
%
%% Syntax
%     WOI = WORDFILT2(I, wk);
%     WOI = WORDFILT2(I, wk, order);
%
%% Inputs
% *|I|* : input image, possibly multispectral.
%
% *|wk|* : kernel defining the neighbourhood and the local integer weights  
%     used in the estimation of the weighted ordered statistics; a possibly  
%     negative weight indicates the position of the centre pixel (origin of 
%     the offset); otherwise, the kernel should be square-shaped with odd  
%     size and its barycenter will be considered as the centre.
%
% *|order|* : optional argument stating the order statistics to be computed in
%     the sorted set of neighbours specified by the nonzero elements in |wk|;
%     it can be either:
%
% * one among the strings |'med'| or |'mod'| when computing the median or 
%         the mode resp. using the Matlab functions with same name,
% * a non null scalar inferior to the sum of weights in |wk| when computing
%         the |order|th element;
%
% note in particular that |order='min'| and |order=1| are equivalent as
%     well as |order='max'| and |order=sum(wk(:))|; however |order='med'|
%     and |order=sum(wk(:))/2| are not equivalent as in the latter case the 
%     value of the median pixel (and not the median value itself like when
%     using |MEDIAN|) is output; default: |order='med'|.
%      
%% Output
% *|WOI|* : output weighted ordered statistics image with same size as the 
%     image |I|; each pixel in |WOI| contains the order statistics in the 
%     neighbourhood of |I| defined by the non-negative values of |wk| around 
%     corresponding pixel in |I|. 
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
% [Soille02]  P. Soille: "On morphological operators based on rank filters",
%       Pattern Recognition, 35:527-535, 2002.
%       <http://www.sciencedirect.com/science/article/pii/S0031320301000474>
%
% [Yaros08]  L. Yaroslavsky: "Local criteria and local adaptive filtering
%       in image processing: a retrospective view", Proc. LNLAIP, 2008.
%       <http://www.eurasip.org/Proceedings/Ext/LNLA2008/papers/cr1005.pdf>
%
%% See also
% Related:
% <matlab:webpub(whichpath('MEDFILT2')) |MEDFILT2|>,
% <matlab:webpub(whichpath('ORDFILT2')) |ORDFILT2|>,
% <WMEDFILT2.html |WMEDFILT2|>.
% Called:
% <WORDFILT2_BASE.html |WORDFILT2_BASE|>.

%% Function implementation
function WOI = wordfilt2(I,wk,varargin)

error(nargchk(2, 4, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

%%
% parsing and checking parameters
if ~isnumeric(I)
    error('wordfilt2:inputerror','matrice is required in input');
end

p = createParser('WORDFILT2');   % create an instance of the inputParser class.
% optional parameters
p.addOptional('order','med', @(x) (isscalar(x) && x>0) ||...
    (ischar(x)&& any(strcmpi(x,{'med','mod'})))); 

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% internal variables and further checking

C = size(I,3);
[x y c] = size(wk);
s = sum(abs(wk(:)));

if ~isequal(round(wk),wk)
    warning('wordfilt2:inputwarning',...
        'the input weighting matrix must have integer values');
end

if isscalar(p.order) && p.order>s
    error('wordfilt2:inputerror',...
        'the input rank must be <= to the sum of the elements in wk');
end

if c~=1 && c~=C
    error('wordfilt2:inputargument', ...
    'input kernel must have 1 channel or the same number as the input');
end

cpos = find(wk<0,1,'first');
if isempty(cpos) && (x~=y || mod(x,2)~=1) 
    error('wordfilt2:inputargument',...
        'input kernel must be square shaped with odd size');
end

%%
% main processing

WOI = wordfilt2_base(I, wk, p.order);

%%
% display

if p.disp
    figure, imagesc(rescale(WOI,0,1)), axis image off;
    if C==1,  colormap gray; end;
    title('weighted order-statistic filtered image');
end

end % end of wordfilt2
