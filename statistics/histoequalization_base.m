%% HISTOEQUALIZATION_BASE - Base function for HISTOEQUALIZATION.
%
%% Syntax
%   x = HISTOEQUALIZATION_BASE(x, y, dir, absval, hinterp);
%   x = HISTOEQUALIZATION_BASE(x, 'gauss', dir, absval, hinterp);
%   x = HISTOEQUALIZATION_BASE(x, 'lin', dir, absval, hinterp);
%
%% See also
% Related:
% <HISTOEQUALIZATION.html |HISTOEQUALIZATION|>,
% Called:
% <matlab:webpub(whichpath('INTERP1')) |INTERP1|>,
% <matlab:webpub(whichpath('SORT')) |SORT|>,
% <matlab:webpub(whichpath('LINSPACE')) |LINSPACE|>.

%% Function implementation
function x = histoequalization_base(x, y, dir, absval, hinterp)

if iscell(x)
    for i=1:length(x)
        if ischar(y),  yy =y;  else yy = y{i};  end
        x{i} = histoequalization_base(x{i}, yy, dir, absval, hinterp);
    end
    return;

elseif ~isempty(dir) && ischar(dir)
    if strcmpi(dir,'col')
        for i=1:size(x,2)
            if ischar(y),  yy =y;  else yy = y(:,i);  end
            x(:,i) = histoequalization_base(x(:,i), yy, dir, absval, hinterp);
        end
        
    elseif strcmpi(dir,'row')
        for i=1:size(x,1)
            if ischar(y),  yy =y;  else yy = y(i,:);  end
            x(i,:) = histoequalization_base(x(i,:), yy, dir, absval, hinterp);
        end
        
    elseif strcmpi(dir,'vec')
        for i=1:size(x,3)
            if ischar(y),  yy =y;  else yy = y(:,:,i);  end
            x(:,:,i) = histoequalization_base(x(:,:,i), yy, dir, absval, hinterp);
        end
    end
    
    return;
end

sx = size(x);
x = x(:);

if ischar(y)
    if strcmpi(y,'gauss')
        y = randn(length(x));
        
    elseif strcmpi(y,'lin') % uniform density
        y = linspace(0,1,length(x));
    
    else
        error('histoequalization_base:inputerror', 'unknown density');
    end
end

y = y(:);

if absval
    s = sign(x);
    x = abs(x);
    y = abs(y);
end

[~,Ix] = sort(x);
vy = sort(y);

nx = length(x);
ny = length(y);

ax = linspace(1,ny,nx)';
ay = (1:ny)';
vx = interp1(ay,vy,ax);
x(Ix) = (1-hinterp)*x(Ix) + hinterp*vx;

if absval,  x = x .* s;  end

x = reshape(x,sx);

end % end of histoequalization_base
