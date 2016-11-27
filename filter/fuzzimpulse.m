%% FUZZIMPULSE - Fuzzy impulse noise reduction.
%
%% Description
% Perform the fuzzy impulse noise reduction in images following the grayscale
% approach of [SWNWK07].
%
%% Syntax
%     F = FUZZIMPULSE(I);
%     [F, N] = FUZZIMPULSE(I, K, L);
%
%% Inputs
% *|I|* : input image, possibly multispectral (then, channels are processed
%     separately).
%
% *|K|* : half-size of the analyzing window used for impulse noise detection
%     (see Section 2.1 of [SWNWK07], Eq.(2)); default: |K=1|.
%
% *|L|* : half-size of the analyzing window used for input image filtering
%     (see Section 3 of [SWNWK07], Eq.(9)); default: |L=1|.
% 
%% Outputs
% *|F|* : output filtered image (impulse noise reduced).
%
% *|N|* : (optional) boolean image set to true for pixels where impulse noise
%     has been detected.
%
%% Remark
% In the case of multispectral inputs, this function applies uncorrelated
% filtering (channel by channel), thus the filtered output may present
% artifact values.
%
%% References
% [SWNWK06]  S. Schulte, V. de Witte, M. Nachtegael, D. van der Weken and
%      E. Kerre: "Fuzzy two-step filter for impulse noise reduction from 
%      color images", IEEE. Trans. on Image Processing, 15(11):3568-3579, 
%      2006.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1709999&tag=1>
%
% [SWNWK07]  S. Schulte, V. de Witte, M. Nachtegael, D. van der Weken and
%      E. Kerre: "Fuzzy random impulse noise reduction method", Fuzzy Sets
%      and Systems, 158(3):270-283, 2007.
%      <http://www.sciencedirect.com/science/article/pii/S0165011406004179>
% 
%% See also
% Related:
% <../../sharpen/html/RAMPSHARP.html |RAMPSHARP|>.
% Called: 
% <FUZZIMPULSE_BASE.html |FUZZIMPULSE_BASE|>. 

%% Function implementation
function [F,varargout] = fuzzimpulse(I, varargin)

error(nargchk(1, 13, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

if ~isnumeric(I)
    error('fuzzimpulse:inputparameter','matrix required in input'); 
end

%% 
% parsing parameters

p = createParser('FUZZIMPULSE');   
p.addOptional('K',1, @(x)isscalar(x) && x>=1 && x==round(x)); 
p.addOptional('L',1, @(x)isscalar(x) && x>=1 && x==round(x)); 

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% main computation

[F,N] = fuzzimpulse_base(I, p.K, p.L);

if nargout==2,  varargout{1} = N;  end; % N computed at no cost

%%
% display

if p.disp
    figure, 
    subplot(1,2,1), imagesc(N), axis image off, title('detected impulse noise'),
    subplot(1,2,2), imagesc(rescale(F)), axis image off, title('denoised output');
    if size(F,3)==1,  colormap gray;  end
end

end % end of fuzzimpulse
