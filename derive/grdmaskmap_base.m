%% GRDMASKMAP_BASE - Base function for GRDMASK.
% 
%% Syntax
%   [Gx,Gy] = GRDMASKMAP_BASE(I, map, method, oaxis)
%
%% See also  
% GRDMASK:
% <GRDMASK.html |GRDMASK|>,
% <GRDMASK_BASE.html |GRDMASK_BASE|>.
% Called:
% <matlab:webpub(whichpath('FSPECIAL')) |FSPECIAL|>.

%% Function implementation
function [Gx,Gy] = grdmaskmap_base(I, map, method, oaxis)

if all(map==true)
    [Gx,Gy] = grdmask_base(I, method, oaxis);
    return
else
    np = sum(map(:));
end

%% 
% dealing with multispectral image

C = size(I,3);

if C>1
    Gx = zeros(np,1,C);    Gy = zeros(np,1,C);
    for ic=1:C
        [Gx(:,:,ic), Gy(:,:,ic)] = grdmaskmap_base(I(:,:,ic), map, method, oaxis);
    end
    return
end

%% 
% main computation

switch method
 
    % define the directional masks
    case {'sobel', 'sob'}
        Mx = -fspecial('sobel'); % [-1 -2 -1; 0 0 0; 1 2 1]
        
    case {'prewitt', 'prew'}
        Mx = -fspecial('prewitt');  % [-1 -1 -1; 0 0 0; 1 1 1]
        
    case {'kirsch', 'kir'}
        Mx = [-5 -5 -5; 3 0 3; 3 3 3];

    case {'robinson', 'rob'}
        Mx = [-1 -1 -1; 1 -2 1; 1 1 1];

    case {'circular', 'circ'}
        Mx = [-0.464 -0.959 -0.464; 0 0 0; 0.464 0.959 0.464];
        
    case {'optimal', 'opt'}
        Mx = [-0.112737 -0.274526 -0.112737; 0 0 0; 0.112737 0.274526 0.112737];
        
    case {'orientation', 'ori'} % optimal filter for orientation
        Mx = [-0.0938 -0.3125 -0.0938; 0 0 0; 0.0938 0.3125 0.0938];
        
    case {'isotropic', 'iso'}
        Mx = [-1 -sqrt(2) -1; 0 0 0; 1 sqrt(2) 1];
        
    case 'roberts'
        M45 = [1 0; 0 -1]; % 45 deg edge responses
        M135 = [0 1; -1  0]; % 135 deg edge responses
        Mx = M45; My = M135;   %#ok
   
    % or compute directly   
    case {'diff', 'difference'} % central difference
        Mx = [0 -1 0; 0 0 0; 0 1 0];

    case {'backward', 'back'}
        Mx = [0 -1 0; 0 1 0; 0 0 0];
        
    case {'forward', 'for'}
        Mx = [0 0 0; 0 -1 0; 0 1 0];
end

I = padarray(I,[1 1],'both', 'replicate');
map = padarray(map,[1 1], false, 'both');
ix = find(map);

[ic,icd] = ixneighbours(I,ix,8);
icd = sortrows([ic,icd],1);
icd = reshape(icd(:,2),[8, np]);
% note that ix=unique(dum) and that np=length(ix)

A = [I(ix)'; I(icd)];

Mx = Mx / sum(abs(Mx(:)));
tmp = Mx(5); Mx(2:5) = Mx(1:4); Mx(1) = tmp;
My = Mx';

Gx = sum(A .* repmat(Mx(:),[1 np]));

Gy = sum(A .* repmat(My(:),[1 np]));

if strcmpi(oaxis,'xy')
    tmp = Gx; Gx = Gy; Gy = -tmp;
end


end % end of grdmaskmap_base

