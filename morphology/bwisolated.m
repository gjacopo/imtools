%% BWISOLATED - Extract isolated pixels from a logical (binary) image.
% 
%% Syntax
%     ISO = BWISOLATED(M);
%
%% Input
% * |M|* : binary (logical) map with size |(X,Y)|.
%
%% Output
% * |ISO|* : map of isolated pixels with size |(2*X-1,2*Y-1)|. 
% 
%% See also 
% Related:
% <matlab:webpub(whichpath('bwmorph')) |BWMORPH|>,
% <matlab:webpub(whichpath('bwthinupsample')) |BWTHINUPSAMPLE|>.
% Called:
% <matlab:webpub(whichpath('bwmorph')) |BWMORPH|>,
% <matlab:webpub(whichpath('conv2')) |CONV2|>.

%% Function implementation
function iso = bwisolated(map)

if ~isempty(ver('images'))
    iso = map & ~bwmorph(map,'clean');
else
    iso = map & conv2(double(e), ones(3,3), 'same')==1;
end

end % end of bwisolated
