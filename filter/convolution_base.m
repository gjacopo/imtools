%% CONVOLUTION_BASE - Base function for CONVOLUTION.
%
%% Syntax
%     y = CONVOLUTION_BASE(x, h, bound);
%
%% See also
% Related:
% <CONVOLUTION.html |CONVOLUTION|>.
% Called:
% <matlab:webpub(whichpath('CONV')) |CONV|>,
% <matlab:webpub(whichpath('CONV2')) |CONV2|>,
% <matlab:webpub(whichpath('FFT')) |FFT|>,
% <matlab:webpub(whichpath('IFFT')) |IFFT|>,
% <matlab:webpub(whichpath('FFT2')) |FFT2|>.

%% Function implementation
function [y,nd] = convolution_base(x, h, bound)

%%
% dealing with multispectral image

if iscell(x)
    y = cell(length(x),1);
    for i=1:length(x)
        [y{i},nd] = convolution_base(x{i}, h, bound);
    end
    return;
end

C = size(x,3);

if C>1
    y = x;
    for ic=1:C
        [y(:,:,ic),nd] = convolution_base(x(:,:,ic), h, bound);
    end
    return
end

%%
% setting the variables

n = size(x);
p = size(h);

% bound = lower(bound);

nd = ndims(x);
if size(x,1)==1 || size(x,2)==1
    nd = 1;
end
if nd==1 
    n = length(x);
    p = length(h);
end

%%
% main computation

switch bound
    
    case 'sym'  % symmetric boundary conditions
        d1 = floor( p/2 );  % padding before
        d2 = p-d1-1;            % padding after
        
        if nd==1      % in 1D 
            x = x(:); h = h(:);
            xx = [ x(d1:-1:1); x; x(end:-1:end-d2+1) ];
            y = conv(xx,h);
            y = y(p:end-p+1);
        elseif nd==2   % in 2D 
            % double symmetry
            xx = x;
            xx = [ xx(d1(1):-1:1,:); ...
                xx; ...
                xx(end:-1:end-d2(1)+1,:) ];
            xx = [ xx(:,d1(2):-1:1),  xx,  xx(:,end:-1:end-d2(2)+1) ];
            
            y = conv2(xx,h);
            y = y( (2*d1(1)+1):(2*d1(1)+n(1)), (2*d1(2)+1):(2*d1(2)+n(2)) );
        end
        
    case 'per'  % periodic boundary conditions
        if p>n
            error('convolution_base:inputerror', ...
                'h filter should be shorter than x.');
        end
        d = floor((p-1)/2);
        if nd==1
            x = x(:); h = h(:);
            h = [ h(d+1:end); zeros(n-p,1); h(1:d) ];
            y = real( ifft( fft(x).*fft(h) ) );
        else
            h = [ h(d(1)+1:end,:); ...
                zeros( n(1)-p(1),p(2) ); ...
                h( 1:d(1),: ) ];
            h = [ h(:,d(2)+1:end),  zeros(n(1),n(2)-p(2)),  h(:,1:d(2)) ];
            y = real( ifft2( fft2(x).*fft2(h) ) );
        end
        
    otherwise
        error('convolution_base:methoderror', ...
            ['method ' bound ' not implemented yet']);
        
end

end % end of convolution_base
