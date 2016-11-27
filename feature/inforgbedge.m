%% INFORGBEDGE - Visualization of edge magnitude and orientation of color images.
%
%% Description
% Provide edge magnitude and orientation of multispectral images based on 
% their derivatives, using multichannel approach (gradient structure
% tensor), and reduced (RGB or LAB) approach (gradient). Used for simple
% visual comparison of both approaches.
%
%% Syntax
%     [E, M, O, rE, rM, rO] = INFORGBEDGE(I, 'disp', true);
%     [E, M, O, rE, rM, rO] = INFORGBEDGE(I, reduce, rho, sig, der, int);
%
%% Inputs
% *|I|* : an input image with size |(X,Y,C)|, where |C=3| when |I| is multichannel.
%   
% *|reduce|* : optional string setting the method used for reducing the input
%      image into a scalar (graylevel) one; it is either |'rgb'| for calling
%      |RGB2GRAY| or |'lab'| for calling |RGB2LAB|; default: |reduce='rgb'|.
%      
% *|rho|* : optional integration scale (see |GTSMOOTH|); default: |rho=0|, 
%      ie. in order to compare both multichannel and reduced approaches, we
%      skip in principle the spatial averaging of the tensor.
%      
% *|sig|* : optional differentiation scale (see |GTSMOOTH|); default: |sig=1|,
%      i.e. Gaussian regularisation is used for estimating the derivatives.
%      
% *|der|* : optional string for pre-smoothing/differentiation (see |GRDSMOOTH|):
%      |'matlab'|, |'vista'|, |'fast'|, |'conv'|, |'fleck'|, |'tap5'|, |'tap7'|,
%      |'sob'|, |'opt'| or |'ana'|; default: |der='matlab'|.
%      
% *|int|* : optional string for post-smoothing (see |GRD2GST|): |'matlab'|,
%      |'conv'| or |'fast'|; default: |int='fast'|.
%
%% Outputs
% *|E, rE|* : edge maps obtained through multichannel and reduced approaches.
%      
% *|M, O|* : magnitude and orientation using the multichannel approach.
%      
% *|rM, rO|* : ibid reducing the image first.
%
% Displays the overlapped gradient magnitudes and edge maps of a color
% image using both multispectral and reduced approaches.
%
%% See also
% Related:
% <EDGECORNER.html |EDGECORNER|>,
% <KOETHEDGE.html |KOETHEDGE|>,
% <CANNYEDGE.html |CANNYEDGE|>.
% Called: 
% <INFORGBEDGE_BASE.html |INFORGBEDGE_BASE|>.

%% Function implementation
function [E, M, O, rE, rM, rO] = inforgbedge(I, varargin)

%% 
% parsing parameters

error(nargchk(1, 19, nargin, 'struct'));
error(nargoutchk(0, 4, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(I)
    error('gstsmooth:inputerror','a matrix is required in input'); 
end

% optional parameters
p = createParser('GSTSMOOTH');   
% principal optional parameters
p.addOptional('reduce', 'rgb', @(x)ischar(x) && any(strcmpi(x,{'rgb','lab'})));
p.addOptional('rho', 0, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addOptional('sigma', 1, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addOptional('der', 'matlab', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','vista','fast','conv','fleck', ...
    'tap5','tap7','sob','opt','ana'}))));
p.addOptional('int', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','conv','fast'})))); % 'ani' not allowed here

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% checking/setting variables

if size(I,3)~=3
    error('inforgbedge:inputerror', 'matrix (m x n x 3) required in input');
end

I = rescale(I,0,1);

if strcmpi(p.reduce,'rgb')
    R = rgb2gray(I);
elseif strcmpi(p.reduce,'lab')
    C = makecform('srgb2lab');
    R = applycform(I,C);
    R = R(:,:,1);
    % R = rgb2lab(I);
end

%%
% main calculations

[E, M, O, rE, rM, rO] = inforgbedge_base(I, R, p.rho, p.sigma, p.der, p.int);

%%
% display

if p.disp
    %  figure, colormap gray
    %  subplot(1,2,1), imagesc(mag), axis image off, title('multichannel')
    %  subplot(1,2,2), imagesc(rmag), axis image off, title('reduced')
    figure,
    subplot(1,2,1), imagesc(rescale(cat(3, M, rM, rM),0,1))
    title('magnitude: R multichannel - GB reduced'), axis image off
    subplot(1,2,2), imagesc(rescale(cat(3, E, rE, rE),0,1))
    title('edge map: R multichannel - GB reduced'), axis image off
end

end % end of inforgbedge
