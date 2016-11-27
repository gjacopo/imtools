function mlf = lifetimeseg(hseg)
% LIFETIMESEG - Pixelwise life time index computed over a stack of labelled
% images (eg., the different levels of a hierarchical segmentation stored in
% increasing order - from fine to coarse - in an array of matrices).
%
%           mlf = lifetimeseg(hseg);
% 
% Input:
%   hseg : an array (X x Y x ns) of matrices sotring a stack of label images, 
%     typically obtained as the output of a hierarchical segmenation, and. 
%     ordered in incresing order, from fine to coarse (ie, each connected
%     region of hseg(:,:,i) belong a (possibly larger) connected region
%     in hseg(:,:,i)).
% 
% Output:
%   mlf : a matrix of size (X x Y) storing the lifetime of every single pixel
%     in the input hierarchical segmentation; 1<=mlf<=ns.
%
% References:
%   [Soille08]  P. Soille: "Constrained connectivity for hierarchical image
%       partitioning and simplification", IEEE Trans. on Pattern Analysis 
%       and Machine Intelligence 30:1132?1145, 2008.
%   [SG09]  P. Soille and J. Grazzini: "Constrained connectivity and 
%       transition regions", Proc. of ISMM, LNCS 5720, pp. 59?69, 2009.
%
% See also 
% Related: 

if ~isnumeric(hseg)
    error('lifetimeseg:inputparameter','a matrix is required in input'); 
end

ns = size(hseg,3);
mlf = zeros(size(hseg(:,:,1)));

if size(hseg,3) == 1
    warning('lifetimeseg:irrelevantentry',...
        'the input should be an array of ns>1 matrices');
    return; % return mlf equal to 1 everywhere
end

tmp = mlf;
for i = 2:ns
    lastmore = hseg(:,:,i) == hseg(:,:,i-1);
    tmp(lastmore) = tmp(lastmore) + 1;
    tmp(~lastmore) = 1;
    increase = lastmore .* tmp>mlf;
    mlf(increase) = mlf(increase) + 1;
end

end