%% REGIONCLEAN_BASE - Base function for REGIONCLEAN.
%
%% Syntax
%   [S, Ck, ColCk] = ...
%         REGIONCLEAN_BASE(L, Ck, ColCk, I, conn, m, features, thresholds);
%
%% See also
% Related:
% <REGIONCLEAN.html |REGIONCLEAN|>,
% <matlab:webpub(whichpath('REGIONPROPS')) |REGIONPROPS|>,
% <IMLABEL.html |IMLABEL|>.
% Called:
% <REGIONADJACENCY_BASE.html |REGIONADJACENCY_BASE|>,
% <matlab:webpub(whichpath('REGIONPROPS')) |REGIONPROPS|>,
% <matlab:webpub(whichpath('BWLABEL')) |BWLABEL|>.

%% Function implementation
function [S, Ck, ColCk] = ...
    regionclean_base(S, Ck, ColCk, I, conn, m, features, thresholds)

%%
% check/set parameters
if ~isequal(length(features),length(thresholds))
    error('regionclean_base:inputerror', ...
        'vectors of feature names and threshold values must have same length' );
elseif ~all(ismember(features,{'Area','Solidity','Extent'}))
    % 'Perimeter','Eccentricity'
    error('regionclean_base:inputerror', ...
        'unknown feature - see REGIONCLEAN and REGIONPROPS helps');
end

[X,Y] = size(S);

%% 
% initialize the labels

% if there is a label of 0 ensure we do not renumber that region by removing
% it from the list of labels to be renumbered
labels = unique(S(:)); 
% labels = sort(S(:)); labels(labels(1:end-1)==labels(2:end)) = []; % inline UNIQUE
labels = setdiff(labels,0);
nlabels = max(labels);
% % cc = imlabel(S,8); % we can extract the connected components
% % ncc = max(cc(:)); 

% if ~isequal(nlabels,length(labels))
%     S = compressregions(S, labels);
%     labels = setdiff(unique(S(:)),0); % update
% end
% %   %----------------------------------------------------------------------
%     function L = compressregions(L, labels)
%         % this function ensures that the labels are set in sequential order
%         % 0 values in the original image are left with the label 0
%         if labels(1)~=0, labels = [0; labels];  end
%         A = cumsum(diff(labels)-1);
%         labels = labels(2:end);
%         A = labels - A;
%         % do the relabeling: compress
%         for l=1:length(labels),    L(L==labels(l)) = A(l);  end
%     end
% %   %----------------------------------------------------------------------

%%
% break regions by assigning a new label to unconnected sub-regions 
% checking that there is only one region for each segment label; if there
% is more than one region they are given unique labels
[S, ischanged] = splitregions(S, conn, labels, nlabels);

%   %----------------------------------------------------------------------
    function [S, ischanged] = splitregions(S, conn, labels, nlabels)
        nl = nlabels;  ischanged = [];
        for l = labels'
            % define the number of connected objects found for each label
            [bl,num] = bwlabel(S==l, conn);  % use given connectedness
            if num > 1  % more than one region with the same label
                Ibl = bl>1; % therefore bl(Ibl)>=2
                S(Ibl) = bl(Ibl) + nl - 1;
                ischanged = [ischanged; l];                            %#ok
                nl = nl + (num - 1); % updated nlabels
            end
        end
        if nl>nlabels,  ischanged = [ischanged; (nlabels+1:nl)'];  end
        % note that we do keep the initial labelling of most of the regions
        % as we add new regions with labels starting after the larger label
        % found in the input image
    end
%   %----------------------------------------------------------------------

if ~isempty(ischanged)  % possibly update
    labels = setdiff(unique(S(:)),0);
    nlabels = max(labels);
end

%%
% update (or not) the centroids of the regions
if isempty(Ck) || ~isempty(ischanged)
    % here, we recompute all centroids for simplicity...
    Ck = nan(nlabels,2);  ColCk = nan(nlabels,3);
    [Ck(labels,:), ColCk(labels,:)] = updatecentroids(S, I, labels);
    % elseif ~isempty(ischanged)
    %  [Ck(ischanged,:), ColCk(ischanged,:)] = updatecentroids(S, I, ischanged);
end

%   %----------------------------------------------------------------------
    function [Ck, ColCk] = updatecentroids(Q, I, k)
        if nargin==3 && ~isempty(k)  % get segments' centroids
            props = regionprops(Q .* ismember(Q,k),'centroid');
        else % perform for all
            props = regionprops(Q,'centroid');
        end
        % note in calling REGIONPROPS: non represented labels in the input 
        % matrix Q are assigned NaN values into the corresponding cells of
        % the calculated feature(s)
        Ck = floor(cat(1, props.Centroid));
        if nargin<3 || isempty(k),  k = 1:size(Ck,1);  end
        Ck = fliplr(Ck(k,:)); % keep only those we are interested in
        % Ck = fliplr(floor(cat(1, props.Centroid))); % keep even NaN values
        cind = Ck(:,1) + ((Ck(:,2)-1)*Y);  % A = isnan(cind);
        cind = [cind cind+X*Y cind+2*X*Y];  % cind(A,:) = [];
        ColCk = nan(size(Ck,1),3); % default output
        if ~isempty(I),
            ColCk = reshape(I(cind), [size(Ck,1) 3]);
            % ColCk(~A,:) = reshape(I(cind), [sum(~A) 3]);
        end
    end
%   %----------------------------------------------------------------------

%%
% get the sparse adjacency matrix for all represented labels
G = regionadjacency_base(S, labels, conn);
% figure, imagesc(label2rgb(S))

%%
% find regions whose features do not respect the given thresholds
R = criteriaregions(S, features, thresholds);

%%
% simple case: get rid of completely isolated regions
isdeleted = find(R);
if isempty(isdeleted)
     % nothing more to do: none of the current regions desobey the given 
     % criterion
     return;
end

%   %----------------------------------------------------------------------
    function C = criteriaregions(S, features, thresholds, k)
        if nargin==4 && ~isempty(k)  % get segment(s)' area(s)
            props = regionprops(S.*ismember(S,k),features);
        else
            props = regionprops(S,features);
        end
        % initialize with the first criterion
        C = cat(1,props.(features{1}));
        % C = C<thresholds(1) & ~isnan(C);
        if nargin<4 || isempty(k),  k = 1:length(C);  end
        C = C(k,:)<thresholds(1) & ~isnan(C(k,:));
        % complete with other possible criteria
        for ip=2:numel(features)
            %if ~isfield(props,features{ip})
            % eval([genvarname(features{ip}) '= cat(1,props.(features{ip}));']);
            A = cat(1,props.(features{ip}));
            C = C | (A(k,:)<thresholds(ip) & ~isnan(A(k,:)));
            % end
        end
    end
%   %----------------------------------------------------------------------

%%
% fill holes by cleaning completely surrounded regions
[S, ischanged] = fillregions(S, G, isdeleted);

%   %----------------------------------------------------------------------
    function [S, ischanged] = fillregions(S, G, isolated)
        ischanged = [];
        % remove isolated regions embodied in another (unique) region
        for l=isolated'
            ind = find(full(G(l,:)));
            if length(ind)==1
                ischanged = [ischanged; ind];                          %#ok
                S(S==l) = ind;
            end
        end
    end
%   %----------------------------------------------------------------------
%     function S = fillholes(S)
%         % remove isolated pixels embodied in another (unique) region
%         isolated = ismember(S,find(area==1 & ~isnan(area)));
%         [ic,icd] = ixneighbours(S,isolated);
%         P = accumarray(ic, S(icd), [], @(x){x});
%         A = cellfun(@(p) length(p)<1,P);
%         P(A) = [];
%         A = find(~A);
%         I = cellfun(@(x) sum(diff([x(end);x(:)],1))==0, P);
%         S(A(I)) = P{I}(1);
%     end
% %   %----------------------------------------------------------------------

%%
% update if any change
if ~isempty(ischanged)
    % % recompress
    % S = compressregions(S, unique(S(:)));
    % update the list of available labels
    labels = setdiff(unique(S(:)),0);
    nlabels = max(labels);
    % update all the representative centroids
    Ck = nan(nlabels,2);  ColCk = nan(nlabels,3);
    [Ck(labels,:), ColCk(labels,:)] = updatecentroids(S, I, labels);
    % update the sparse adjacency matrix
    G = regionadjacency_base(S, labels, conn);
    % update the area vector
    R = criteriaregions(S, features, thresholds);
    isdeleted = find(R);
end

s = floor(sqrt(X*Y/length(labels))); % kind of arbitrary...

%%
% assign those regions the label of the adjacent region whose centroid is
% closest to the centroid of the current region
for k = isdeleted' 
    % k is a label index of a region we will have to merge to another one
    r = true; % as said, we will merge it for sure    
    % find regions adjacent to k
    ind = find(full(G(k,:)));
    
    %%
    % keep merging with the closest element in the adjacency matrix until
    % we obtain an area >= tharea, or we run out of regions to merge
    while ~isempty(ind) && r
        %%
        % compute the distances between the centroid of the current
        % region and the centroids of all its neighbours
        if length(ind)>1
            ds = calclabspace(Ck(repmat(k,size(ind)),:), Ck(ind,:), ...
                ColCk(repmat(k,size(ind)),:), ColCk(ind,:), s, m);
            [~,i] = min(ds);
        else
            i = 1; % one neighbour only
        end
        %%
        % merge both regions k and ind(i) and store the new merged region
        % in k; also update the connecting graph
        [S, G] = updateregions(S, G, k, ind(i));
        %%
        % update regions' centroids for the merged region k only
        [Ck(k,:), ColCk(k,:)] = updatecentroids(S, I, k);
        %%
        % update the criteria on the merged region k only
        r = criteriaregions(S, features, thresholds, k);
        ind = find(full(G(k,:))); % the adjacency matrix has changed
    end
    
end

%%
% note that this is an arbitrary merging, depending in particular of the
% order the different regions are merged to others...

%   %----------------------------------------------------------------------
    function [S, G] = updateregions(S, G, s1, s2)
        % function to merge segment s2 into s1: the segmentation image, the
        % adjacency matrix and the area vector are updated.
        
        % update the label image: relabel s2 with s1 in the segment image
        S(S==s2) = s1;
        
        % update the connecting graph: s1 inherits the adjacancy matrix
        % entries of s2
        G(s1,:) = G(s1,:) | G(s2,:);
        G(:,s1) = G(:,s1) | G(:,s2);
        G(s1,s1) = 0;  % ensure s1 is not connected to itself
        % disconnect s2 from the adjacency matrix
        G(s2,:) = 0;
        G(:,s2) = 0;
    end
%   %----------------------------------------------------------------------
    function Ds = calclabspace(p0, p1, Lab0, Lab1, s, m)
        distL2 = @(v0, v1) sqrt(sum(abs(v1-v0).^2, 2));
        
        dxy =  distL2(p0, p1);
        dlab = distL2(Lab0, Lab1);
        % Ds is the sum of the lab distance and the xy plane distance
        % normalized by the grid interval S.
        if isempty(dlab),  dlab=0;  end
        if isempty(dxy),   dxy=0;  end
        Ds = dlab + m * dxy / s;
        % note: the greater the value of m, the more spatial proximity is
        % emphasized and the more compact the cluster.
    end
%   %----------------------------------------------------------------------

% % as some regions will have been absorbed into others and no longer exist
% % we now renumber the regions so that they sequentially increase from 1
% labels = unique(S(:));  % sorted list of unique labels
% S = compressregions(S, labels);
% % Ck and ColCk unchanged

end % end of regionclean_base
