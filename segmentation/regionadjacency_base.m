%% REGIONADJACENCY_BASE - Base function for REGIONADJACENCY.
%
%% Syntax
%      G = REGIONADJACENCY_BASE(L, labels, conn);
%
%% Remark
% Note that the Image Processing toolbox is not mandatory; in the case it
% is available, it will however systematically be used.
%
%% See also  
% Related:
% <REGIONADJACENCY.html |REGIONADJACENCY|>.
% Called:
% <matlab:web(whichpath('SPARSE')) |SPARSE|>,
% <matlab:web(whichpath('SORT')) |SORT|>,
% <matlab:web(whichpath('UNIQUE')) |UNIQUE|>,
% <matlab:web(whichpath('IMDILATE')) |IMDILATE|>,
% <matlab:web(whichpath('DIFF')) |DIFF|>,
% <matlab:web(whichpath('SETDIFF')) |SETDIFF|>.

%% Function implementation
function G = regionadjacency_base(L, labels, conn)

labels = labels(:);

% if isempty(labels),  labels = setdiff(unique(L(:)),0);  end
nlabels = max([labels; unique(L(:))]);

if ~isempty(ver('images'))
    
    I = [];
    
    if conn==4,      SE = strel('diamond',1);
    elseif conn==8,  SE = strel('square',3);
    end
    
    for l = labels'
        % dilate each labeled region l and use it as a mask on the input 
        % label image: this extracts the original region plus a one pixel
        % wide section of any adjacent regions
        r = L(imdilate(L==l, SE));
        % find the list of (unique) labels connected to this region, ensuring 
        % the label l itself and 0 are not included in this list
        a = setdiff(unique(r), [l 0]);
        % form the list of labels that are adjacent to the current region l
        I = [I; repmat(l,[length(a) 1]) a(:)];                            %#ok
    end
    
else
    
    % general 4-connectivity: compare neighbours to the left and to the bottom
    L2 = L(2:end,:); L1 = L(1:end-1,:); % matrices dimension (X-1,Y)
    ind = find(L2~=L1); % find(diff(L,1,1)~=0);
    I = [L2(ind) L1(ind)];
    L2 = L(:,2:end);  L1 = L(:,1:end-1); % matrices dimension (X,Y-1)
    ind = find(L2~=L1); % find(diff(L,1,2)~=0);
    I = [I; L2(ind) L1(ind)];
    
    if conn==8
        % 8-connectivity: compare to the diagonal neighbours
        L2 = L(2:end,2:end);  L1 = L(1:end-1,1:end-1); % matrices dimension (X,Y-1)
        ind = find(L2~=L1);
        I = [I; L2(ind) L1(ind)];
        L2 = L(2:end,1:end-1);  L1 = L(1:end-1,2:end); % matrices dimension (X,Y-1)
        ind = find(L2~=L1);
        I = [I; L2(ind) L1(ind)];
    end
    
    % regions with a label of 0 are not considered.
    I(I(:,1)==0 | I(:,2)==0,:) = [];
end

% get rid of duplicates
I = unique(sort(I,2),'rows');
% form the sparse adjacency matrix
G = sparse(I(:,1), I(:,2), ones(length(I),1), nlabels, nlabels);
G = G | G';

end % end of regionadjacency_base