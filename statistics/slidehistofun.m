%% SLIDEHISTOFUN - Local histogram-based image statistics.
%
%% Description
% Function for computing local statistical features based on local histograms
% computed over square sliding windows.
%
%% Syntax
%    O = SLIDEHISTOFUN(I);
%    O = SLIDEHISTOFUN(I, func, win, bin, method);
%
%% Inputs
% *|I|* : input image with size |(X,Y,C)| (|C>1| for multispectral images)
%     over which local statistics are computed.
%
% *|func|* : function handle or cell of function handles to apply on the local
%     histograms.
%
% *|win|* : (optional) size of the local (square) sliding window over which
%     statistics are computed; default |win=3|.
%
% *|bin|* : (optional) number of bins used in the histogram representation
%     as resolution reduction; default: |bin| is the max number of values
%     available in the input image. 
%
% *|method|* : (optional) string indicating the method choose for computing
%     the local histograms in sliding windows; it is either:
%
% * |'int'| for implementing the integral histogram approach of [Pori05]
%          (see function |HISTOINTEGRAL1D_BASE|),
% * |'joint'| for implementing the joint integral histogram approach of 
%          [ZLLG10] (see function |HISTOINTEGRALJOINT_BASE|),
% * |'dist'| for implementing the distributive histogram approach of [SDH08];
%     
% note that in the two first mentioned cases, the integral histograms need 
%     to be computed and stored prior to the function application; default: 
%     |method='dist'|;
% 
%% Output
% *|O|* : an image of size |(X,Y)| when a single function is computed over 
%     the local histograms (|func| passed as a function handle, see above),
%     or a cell of such images when several functions are estimated (|func|
%     passed as a cell of functions handles); |O(x,y)| (or |O{i}(x,y)|) is
%     the output of the function func (or |func{i}|) of the histogram computed
%     on the local window centered at |(x,y)|. 
% 
%% References
% [Pori05]  F. Porikli: "Integral histogram: a fast way to extract histogram
%      features", Proc. IEEE CVPR, pp. 829-836, 2005.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1467353>
%
% [ZLLG10]  K. Zhang, G. Lafruit, R. Lauwereins and L. Van Gool: "Joint 
%      integral histograms and its application in stereo matching", Proc. 
%      IEEE ICIP, pp. 817-820, 2010.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5653410>
%
% [SDH08]  M. Sizintsev, K.G. Derpanis, A. Hogue: "Histogram-based search:
%      a comparative study", Proc. IEEE CVPR, pp. 1-8, 2008.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4587654>
%
%% See also
% Related:
% <matlab:webpub(whichpath('BLKPROC')) |BLKPROC|>,
% <matlab:webpub(whichpath('COLFILT')) |COLFILT|>.
% Called:
% <SLIDEHISTOFUN_BASE.html |SLIDEHISTOFUN_BASE|>.

%% Function implementation
function O = slidehistofun(I, varargin)

%%
% parsing parameters

narginchk(2, 13);
nargoutchk(1, 1);

if ~isnumeric(I)
    error('slidehistofun:inputparameter','vector or matrix required in input'); 
end

p = createParser('SLIDEHISTOFUN');   
p.addRequired('func', @(x)iscell(x) || ...
    strcmpi(class(x),'function_handle'));  % has to be entered!!!
p.addOptional('win', 3, @(x)isscalar(x) && x>=1);  
p.addOptional('bin', [], @(x)isempty(x) || ...
    (isscalar(x) && x>=1 && round(x)==x));
p.addOptional('method', 'dist', @(x)ischar(x) && ...
    any(strcmpi(x,{'dist','int','joint'})));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% checking variables

if isempty(p.bin),  p.bin = max(I(:)) + 1;  end

if mod(p.win,2)==0,  p.win = p.win + 1;  end

%%
% main computation

O = slidehistofun_base(I, p.method, p.func, p.win, p.bin);

%%
% display

if p.disp
    if iscell(O)
        N = numel(O);
        figure; ndisp = 3; mdisp = ceil(N/ndisp);
        for i=1:N
            subplot(mdisp, ndisp, i),  imagesc(rescale(O{i})), axis image off;
            title(func2str(p.func{1}));
        end
    else
        figure, imagesc(rescale(O)), axis image off, title(func2str(p.func));
    end
    
end

end % end of slidehistofun
