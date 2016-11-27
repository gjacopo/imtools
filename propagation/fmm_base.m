%% FMM_BASE - Base function for FMM.
%
%% Description
%
%% Syntax
%  [D, P] = FMM_BASE(W, m_seed, m_path, start_pts, end_pts, order, L, niter);
%
%% See also
% Related:
% <FMM.html |FMM|>,
% <IM2POTENTIAL.html |IM2POTENTIAL|>,
% <POTENTIAL2FRONT.html |POTENTIAL2FRONT|>,
% <../../graph/html/DIJKSTRA.html |DIJKSTRA|>.
% Called: 
% <FMMISOPROPAGATION.html |FMMISOPROPAGATION|>,
% <FMMANISOPROPAGATION.html |FMMANISOPROPAGATION|>,
% <FMMPATH.html |FMMPATH|>.

%% Function implementation
function [D, P] = ...
    fmm_base(W, m_seed, m_path, start_pts, end_pts, order, L, niter)

%%
% check the input parameters

% we still allow variable nargin for this function
if nargin<=7,    niter = Inf;  
    if nargin<=6,  L = [];  end
end

n = size(start_pts,1);
m = size(end_pts,1);

if islogical(m_path) && ~m_path,  P = [];  end

%% 
% calculation

d = nb_dims(W);

if d==3 || d==4 % vector field or tensor field: anisotropic FMM
    
    if any(strcmpi(m_seed,{'sing','allpairs'}))
        D = fmmanisopropagation(W, start_pts, L);
        if ischar(m_path)
            P = fmmpath(D, start_pts, end_pts, m_path);
        end
        
    elseif strcmpi(m_seed,'mult') % distance from multiple sources
        D = cell(n,1);  
        if ischar(m_path),  P = cell(n,m);  end
        for i=1:n
            D{i} = fmmanisopropagation(W, start_pts(:,i), L);
            if ischar(m_path)
                P(i,:) = fmmpath(D{i}, start_pts(:,i), end_pts, m_path);  %#ok
            end
        end
        
    end
    
elseif d==2 % scalar field: isotropic FMM
    
    if any(strcmpi(m_seed,{'sing','allpairs'}))
        D = fmmisopropagation(W, start_pts, end_pts, order, L, niter);
        if ischar(m_path)
            P = fmmpath(D, start_pts, end_pts, m_path);
        end
        
        
    elseif strcmpi(m_seed,'mult') % distance from multiple sources
        D = cell(n,1);  
        if ischar(m_path),  P = cell(n,m);  end
        for i=1:n
            D{i} = fmmisopropagation(W, start_pts, order, L, niter);
            if ischar(m_path)
                P(i,:) = fmmpath(D{i}, start_pts(:,i), end_pts, m_path);  %#ok
            end
        end
        
    end
    
else
    error('fmm_base:inputerror', ...
        'input cost function should be a vector of a tensor field');
end

end % end of fmm_base



