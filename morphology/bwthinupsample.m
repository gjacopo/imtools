%% BWTHINUPSAMPLE - Upsampled morphological thinning of a binary image.
% 
%% Description
% Perform the morphological thinning of an upsampled binary image. Enable to 
% deal with the indetermination of the thinning algorithm for binary structures
% of even width (in number of pixels).
%
%% Syntax
%     BT = BWTHINUPSAMPLE(M);
%
%% Input
% *|M|* : binary (logical) map with size |(X,Y)|.
%
%% Output
% *|BT|* : thinned map with size |(2*X-1,2*Y-1)|. 
%
%% See also
% Related:
% <BWISOLATED.html |BWISOLATED|>,
% <matlab:web(whichpath('BWMORPH')) |BWMORPH|>,

%% Function implementation
function BT = bwthinupsample(M)

if isempty(ver('images'))
    error('bwthinupsample:errortoolbox', 'Image Processing toolbox required');
end

[X,Y] = size(M);

M2 = zeros(2*X-1,2*Y-1);
M2(1:2:end,1:2:end) = M;

M2(2:2:end,1:2:end) = M(1:end-1,:) & M(2:end,:);
M2(1:2:end,2:2:end) = M(:,1:end-1) & M(:,2:end);

M2(2:2:end,2:2:end) = (M(1:end-1,1:end-1) & M(2:end,2:end)) | ...
    (M(2:end,1:end-1) & M(1:end-1,2:end));

% thin objects to lines: remove pixels so that an object without holes shrinks 
% to a minimally connected stroke, and an object with holes shrinks to a
% connected ring halfway between each hole and the outer boundary.
BT = bwmorph(M2,'thin',Inf);

end % end of bwthinupsample
