%% GRDMASK_BASE - Base function for GRDMASK.
% 
%% Syntax
%   [Gx,Gy] = GRDMASK_BASE(I, method, oaxis)
%
%% See also  
% Related:
% <GRDMASK.html |GRDMASK|>,
% <GRDMASKMAP_BASE.html |GRDMASKMAP_BASE|>.
% Called:
% <matlab:webpub(whichpath('DERGRADIVATIVES')) |GRAD|>,
% <matlab:webpub(whichpath('FSPECIAL')) |FSPECIAL|>,
% <matlab:webpub(whichpath('IMFILTER')) |IMFILTER|>,
% <matlab:webpub(whichpath('GRADIENT')) |GRADIENT|>.

%% Function implementation
%--------------------------------------------------------------------------
function [Gx,Gy] = grdmask_base(I, method, oaxis)

%%
% dealing with multispectral image

C = size(I,3);

if C>1
    Gx = zeros(size(I));    Gy = zeros(size(I));
    for ic=1:C
        [Gx(:,:,ic), Gy(:,:,ic)] = grdmask_base(I(:,:,ic), method, oaxis);
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
    case {'matlab', 'diff', 'difference'} % central difference
        [Gy,Gx] = gradient(I);

    case {'backward', 'back'}
        [Gx,Gy] = grad(I, 1, 'sym');
        Gx = Gx/2; Gy = Gy/2;

    case {'forward', 'for'}
        [Gx,Gy] = grad(I, 2, 'sym');
        Gx = Gx/2; Gy = Gy/2;

    case {'derivative5', 'tap5'}
        [Gy,Gx] = derivative5(I, 'x', 'y');

    case {'derivative7', 'tap7'}
        [Gy,Gx] = derivative7(I, 'x', 'y');

end

if ~any(strcmpi(method,{'matlab','diff','difference', 'backward','back', ...
        'forward','for', 'derivative5','tap5', 'derivative7','tap7'}))
    % define the horizontal mask
    My = Mx';
    
    % filter to get the derivatives
    Gx = imfilter(I,Mx/sum(sum(abs(Mx))),'replicate');  % vertical
    Gy = imfilter(I,My/sum(sum(abs(My))),'replicate'); % horizontal
end

if strcmpi(oaxis,'xy')
    tmp = Gx; Gx = Gy; Gy = -tmp;
end


end % end of grdmask_base


%% Subfunction

%%
% |GRAD| - gradient, forward and backward differences.
%
%    [fx,fy,fz] = grad(M, order, bound);
%    bound : 'per' or 'sym'
%    order : 1 (backward differences) or 2 (forward differences).
%
% Assumes that the function is evenly sampled with sampling step 1.
% Note: the grad operator is *minus* the transpose of the div operator.
%--------------------------------------------------------------------------
function [fx,fy,fz] = grad(M, order, bound)

% retrieve number of dimensions
nbdims = 2;
if size(M,1)==1 || size(M,2)==1
    nbdims = 1;
end
if size(M,1)>1 && size(M,2)>1 && size(M,3)>1
    nbdims = 3;
end


if strcmp(bound, 'sym')    
    if order==1
        fx = M([2:end end],:,:)-M;
    else
        fx = ( M([2:end end],:,:)-M([1 1:end-1],:,:) )/2;
        % boundary
        fx(1,:,:) = M(2,:,:)-M(1,:,:);
        fx(end,:,:) = M(end,:,:)-M(end-1,:,:);
    end
    if nbdims>=2
        if order==1
            fy = M(:,[2:end end],:)-M;
        else
            fy = ( M(:,[2:end end],:)-M(:,[1 1:end-1],:) )/2;
            % boundary
            fy(:,1,:) = M(:,2,:)-M(:,1,:);
            fy(:,end,:) = M(:,end,:)-M(:,end-1,:);
        end
    end
    if nbdims>=3
        if order==1
            fz = M(:,:,[2:end end])-M;
        else
            fz = ( M(:,:,[2:end end])-M(:,:,[1 1:end-1]) )/2;
            % boundary
            fz(:,:,1) = M(:,:,2)-M(:,:,1);
            fz(:,:,end) = M(:,:,end)-M(:,:,end-1);
        end
    end
else
    if order==1
        fx = M([2:end 1],:,:)-M;
    else
        fx = ( M([2:end 1],:,:)-M([end 1:end-1],:,:) )/2;
    end
    if nbdims>=2
        if order==1
            fy = M(:,[2:end 1],:)-M;
        else
            fy = ( M(:,[2:end 1],:)-M(:,[end 1:end-1],:) )/2;
        end
    end
    if nbdims>=3
        if order==1
            fz = M(:,:,[2:end 1])-M;
        else
            fz = ( M(:,:,[2:end 1])-M(:,:,[end 1:end-1]) )/2;
        end
    end
end

if nargout==1
    if nbdims==2
        fx = cat(3,fx,fy);
    elseif nbdims==3
        fx = cat(4,fx,fy,fz);
    end
end
end
