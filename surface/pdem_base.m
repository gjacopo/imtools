 function dem = pdem_base(I, seed, method, N, a, ...
     rho, sigma, der, int, samp, eign)
% PDEM_BASE - Base function for PDEM. 
% 
%        dem = pdem_base(I, seed, method, N, a, ...
%                        rho, sigma, der, int, samp, eign);
%
%
% See also PDEM, GEODESICFILT_BASE.
% Calls IM2POTENTIAL, POTENTIAL2FRONT, GSTSMOOTH_BASE, GSTFEATURE_BASE, GSTDECOMP.

if isempty(a),        a = [1 2];
elseif length(a)==1,  a = [a a+eps];
% elseif length(a)==2 &&  a(1)>a(2),  a = a([2 1]); 
end

[n,m,c] = size(I);                                                     %#ok

%% Main computation

% prelimnary estimation of the structure tensor
% if any(strncmp(p.method, {'tensor','hybrid'},6))
%     h = ceil(3*max(p.sig,p.rho));
%     A = padarray(I,[h h],'replicate','both');
%     T = gstsmooth(A, p.rho, p.sig, 'der','pey', 'sm','pey','samp',1); 
%     L = gstfeature(T(:,:,1,1), T(:,:,2,2), T(:,:,1,2),'norm','eign',p.eign);
%     T = T(h+1:n+h, h+1:h+m,:,:);
%     L = L(h+1:n+h, h+1:h+m,:,:);
%     L = rescale(L,0,1); %L = rescale(L,0,1-eps);
% end

%h = ceil(3*max(sigma,rho));
h=0;
%A = padarray(I,[h h],'replicate','both');
A=I;

% N = rescale(N,0,1);

if strcmp(method,'hybrid0')
    % special case, not implemented in IM2POTENTIAL_BASE
    T = gstsmooth_base(A, rho, sigma, der, int, samp, ...
        [], false, false, 8, .4); % default choices
    L = gstfeature_base(T(:,:,1,1), T(:,:,2,2), T(:,:,1,2), ...
        'norm', eign, [], []);
    
    figure, plot_tensor_field(T,rescale(A));
    
    % L = rescale(N,0,1); % L = rescale(N,0,1-eps);
    [l1, l2, e1, e2] = gstdecomp(T);                                   %#ok
    l1 = 1 ./  (1+N).^a(1);
    % l2 = (1+max(L(:))-L) ./ (1+N).^p.a(1);
    l2 = (1-L/max(L(:))) ./ (1+N).^a(2);
    T = gstdecomp(l1,l2,e1,e2);
    
    %     % special case, not implemented in IM2POTENTIAL_BASE
    %     T = gstsmooth_base(A, rho, sigma, der, int, samp, ...
    %         [], false, false, 8, .4); % default choices
    %     L = gstfeature_base(T(:,:,1,1), T(:,:,2,2), T(:,:,1,2), ...
    %         'norm', eign, [], []);
    %     figure, imagesc(L), colormap gray
    %     % L = rescale(N,0,1); % L = rescale(N,0,1-eps);
    %     [l1, l2, e1, e2] = gstdecomp(T);
    %     l1 = (1-L/max(L(:))) ./  (1+N).^a(2);
    %     % l2 = (1+max(L(:))-L) ./ (1+N).^p.a(1);
    %     l2 = 1 ./ (1+N).^a(1);
    %     T = gstdecomp(l1,l2,e1,e2);
    
elseif strcmp(method,'hybrid2')
    % compute pdem using the isotropic tensor field of the image
    % orthogonal vector
    T = gstsmooth_base(A, rho, sigma, der, int, samp, ...
        [], false, false, 8, .4); % default choices
    V = cat(3, -T(:,:,2), T(:,:,1));
    % new tensor
    T = gstdecomp(ones(n,m), ones(n,m), T, V );
    
elseif strcmp(method,'hybrid3')
    % compute pdem using the unitary tensor: equivalent to euclidean
    % isotropic
    T = zeros(n,m,2,2);
    % unitary eigenvector
    T(:,:,2)=1;             T(:,:,1)=0;
    % orthogonal vector
    V = cat(3, -T(:,:,2), T(:,:,1));
    % new tensor
    T = gstdecomp(ones(n,m), ones(n,m), T, V );
    
else  % all other cases implemented in IM2POTENTIAL_BASE
    
    switch method
        case 'intens',  method = 'pixinv'; 
            % compute the classical PDEM
                        
        case 'hybrid1',  method = 'gstninv'; 
            % compute pdem using the unitary tensor weightened by the input
            % image: equivalent to classical geodesic propagation on the
            % input image using the option method='intensity'
            
        case 'hybrid4',  method = 'gstorth';
            %compute pdem using the unitary tensor: equivalent to euclidean
            % isotropic
        
        case 'hybrid5',  method = 'gst';
           
    end
    
    T = im2potential(A, method, a, rho, sigma, der, int, samp, eign);
    
end

% compute the dem
dem = potential2front(T, seed);
dem = dem(h+1:n+h, h+1:h+m,:,:);

end
