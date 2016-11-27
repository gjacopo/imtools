%% TENSANIFILT_BASE - Base function for TENSANIFILT. 
%
%% Syntax
%     [F, S, T] = TENSANIFILT_BASE(I,varargin);
%
%% See also
% <TENSANIFILT.html |TENSANIFILT|>,
% <CONVOLUTION_BASE.html |CONVOLUTION_BASE|>,
% <MDLFILT_BASE.html |MDLFILT_BASE|>,
% <ADAPTIVEFILT_BASE.html |ADAPTIVEFILT_BASE|>,
% <GEODESICFILT_BASE.html |GEODESICFILT_BASE|>.
% Called:
% <TENSCALEDIRFILT_BASE.html |TENSCALEDIRFILT_BASE|>,
% <../../derive/html/HESSMOOTH_BASE.html |HESSMOOTH_BASE|>,
% <../../derive/html/GSTSMOOTH_BASE.html |GSTSMOOTH_BASE|>,
% <../../derive/html/GSTDECOMP.html |GSTDECOMP|>.

%% Function implementation
function [F, S, T] = tensanifilt_base(I, method, rho, sigma, der, int, ...
    samp, a, c, alpha, beta, eps, p1, p2)

%%
% setting internal variables
C = size(I,3);

% initialize the output
F = I;

if any(strcmp(method,{'tsc','tschumperle'})) && p1 < p2
    tmp = p1;    p1 = p2;    p2 = tmp;
end


% [dl,depth,best] = mdlfilt(I,scales,'int',int,'lam',lam);
% S = scales(best);
S=[];

if any(strcmp(method,{'sol','sole'}))
    % compute the normalized Hessian tensor
    T = hessmooth_base(I, rho, sigma, der, int, 1, [], false, true, 8, .4);
else
    % compute the Gradient structure tensor
    T = gstsmooth_base(I, rho, sigma, der, int, samp, [], false, false, 8, .4);
end


%if any(strcmp(method,{'wei','weickert','kim','kimmel',...
%                        'mid','middendorf','sol','sole'}))

if ~strcmp(method,'gst')
    [l1,l2,e1,e2] = gstdecomp(T);   
    
    % modify the eigenvalues
    switch method
        case {'wei','weickert','kim','kimmel'}
            % approaches from [Weick97] and [KMS00]
            nu = (l1 -l2) .* (l1 - l2); % coherence measure
            ir = nu<eps;
            l1 = a * ones(size(nu));
            l2(~ir) = l2(~ir) + (1-a) .* exp(-c ./ nu(~ir));
            if any(strcmp(method,{'kim','kimmel'}))
                l1 = - 1 ./l2;
            end
             % eigenvectors
             tmp = e1;  e1 = e2;  e2 = tmp;
            
        case {'mid','middendorf'} % approach from [MN02]
            % eigenvalues
            tmp = 1 ./ (sqrt(l1)+eps);
            l1 = 1 ./ (sqrt(l2)+eps);
            l2 = tmp;
            % eigenvectors unchanged
           
        case {'sol','sole'} % approach from [SLS]
            k1 = max(abs(l1),abs(l2));
            k2 = l1 + l2 - k1; % min(abs(l1),abs(l2));
            nu = (k1 -k2) ./ (k1+k2);
            nur = zeros(size(l1));
            nuv = nur;
            ir = k1<0;
            nur(ir) = nu(ir);  nuv(ir) = 0;
            ir = k1>0;
            nur(ir) = 0;  nuv(ir) = nu(ir);
            
            l1 = eps;
            l2 = alpha * nur + beta * nuv;
            
        case {'tsc','tschumperle'}
             tmp = l1;
             l1 = power(1+l1+l2, -p1);
             l2 = power(1+l1+l2, -p2);
             l1 = tmp;
             
        case 'gstort'
             tmp = e1;  e1 = e2;  e2 = tmp;            
            
    end
    
    % recompose the structure tensor with modified eigenvalues and/or
    % eigenvectors
    T = gstdecomp(l1,l2,e1,e2);

end

%%
% perform adaptive filtering
for c=1:C
    F(:,:,c) = tenscaledirfilt_base(I(:,:,c), T, S, 8, 3, 5, 12, 4, 0.5, false);
   %  F(:,:,c) = perform_directional_filtering(I(:,:,c),T);
end

end % end of 
