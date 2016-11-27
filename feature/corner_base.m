%% CORNER_BASE - Base function for CORNER.
% 
%% Syntax
%     [cornermap, ptcorner] = ...
%          CORNER_BASE(I, meth, thr, kap, rad, nnmax, gap, thang, sigma, rho, reduce);
%
%% See also
% Related: 
% <CORNER.html |CORNER|>,
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>.
% Called:
% <HARRISCORNER.html |HARRISCORNER|>,
% <SUSANCORNER_BASE.html |SUSANCORNER_BASE|>,
% <FASTCPDA_BASE.html |FASTCPDA_BASE|>,
% <FASTCORNER_BASE.html |FASTCORNER_BASE|>.

%% Function implementation
function [cornermap, ptcorner] = ...
    corner_base(I, method, thres, kappa, radius, nonmax, gap, thang, sigma, rho, reduce )

%%
% dealing with multispectral images

[X,Y,C] = size(I);

if C>1
    ptcorner = cell(C,1);
    cornermap = false(size(I));
    if isempty(sigma),  sigma=zeros(1,1,3);  % dummy variable
    elseif isscalar(sigma), sigma=repmat(sigma,[1 1 3]);   end
    if isempty(rho),  rho=zeros(1,1,3); % dummy variable
    elseif isscalar(rho), rho=repmat(rho,[1 1 3]);   end;  
    for ic=1:C
        [cornermap(:,:,ic), ptcorner{ic}] = ...
            corner_base(I(:,:,ic),method, thres, kappa, ...
            radius, nonmax, gap, thang, sigma(:,:,ic), rho(:,:,ic));
    end
    if reduce
        for c=2:C
            cornermap(:,:,1) = cornermap(:,:,1) | cornermap(:,:,c);
            ptcorner{1} = [ptcorner{1}; ptcorner{c}];
        end
        cornermap = cornermap(:,:,1);
        ptcorner = ptcorner{1};
    end
    return
end

%%
% calling appropriate method

switch method
    case {'harris','noble','forster'} % Harris
        % use Kovesi implementation if desired
        % [cim, r, c] = harris(I, sigma, 1000, radius);
        % ptcorner{1} = [r,c];
        ptcorner = harriscorner_base(I, sigma, rho, kappa, thres, radius);
        
    case 'susan'
        ptcorner = susancorner_base(I, 'c', thres, 'flat', false, true);
        
        %     case 'cur'
        %         ptcorner = curvecorner_base(I);
        
    case {'fast9','fast10','fast11','fast12'}
        nn = method(5:end);
        ptcorner = fastcorner_base(I, nn, nonmax, thres);
        class(ptcorner)
        
    case 'cpda'
        V = cannyedge_base(I, sigma, 'matlab', [], false);
        ptcorner = fastcpda_base(V, gap, thang, false);
        
        % case 'comp'
        % case 'cong'
    otherwise
        error('corner:errorinput', ...
            ['method ' method ' not implemented yet'])
end

ptcorner = ptcorner{1};

cornermap = false(X,Y);
%   cornermap(sub2ind([X,Y],ptcorner(:,1),ptcorner(:,2))) = true;
cornermap(ptcorner(:,1) + (ptcorner(:,2)-1)*X) = true;

end % end of corner_base

