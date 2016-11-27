%% REGIONADJACENCY - Region adjacency matrix of a labelled image.
%
%% Description
% Compute a (sparse) adjacency matrix for a labelled image of segmentation.
%
%% Syntax
%         G = REGIONADJACENCY(L);
%         G = REGIONADJACENCY(L, labels, conn);
%
%% Inputs
% *|L|* : a segmentation image, where each pixel is assigned the (integer) 
%     label of the region it belongs to.  
%
% *|labels|* : (optional) labels of interest: the graph will be output for
%     those labels only; default: |labels=[]|, ie. all (non null) labels
%     present in the input image will be represented.
%
% *|conn|* : (optional) 4 or 8 connectivity; default: |conn=8|.
%
%% Outputs
% *|G|* : a (sparse) adjacency matrix indicating (through non zero entries) 
%     which regions are adjacent to each other, ie. pixels belonging to
%     those regions share boundaries.
%
%% Remarks
% * Note that regions with a label of 0 are not considered.  
% 
% * Two naive naive/intuitive implementations are provided, one of them using
% the Image Processing toolbox if available (see |REGIONADJACENCY_BASE|).
%
%% See also  
% Related:
% <matlab:web(whichpath('LABELMATRIX')) |LABELMATRIX|>,
% <matlab:web(whichpath('BWLABEL')) |BWLABEL|>,
% <matlab:web(whichpath('BWCONNCOMP')) |BWCONNCOMP|>.
% Called:
% <REGIONADJACENCY_BASE.html |REGIONADJACENCY_BASE|>.

%% Function implementation
function G = regionadjacency(L, varargin)

%%
% parsing parameters

error(nargchk(1, 3, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(L)
    error('regionadjacency:inputerror','image of integer labels required in entry'); 
end

p = createParser('REGIONADJACENCY');   % create an instance of the inputParser class.
p.addOptional('labels', [], @(x)isempty(x) || (isnumeric(x) && all(x(:)>=0)));
p.addOptional('conn', 8, @(x)isscalar(x) && (x==4 || x==8));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p); 

%%
% setting variables

if isempty(p.labels)
    % identify the unique labels in the image, excluding 0 as a label
    p.labels = setdiff(unique(L(:)),0);  
elseif max(p.labels(:))>max(L(:))
    warning('regionadjacency:inputwarning', ...
        'labels not represented in the input segmentation are ignored');
    p.labels(p.labels(:)>max(L(:))) = [];
end


%% 
% main computation

G = regionadjacency_base(L, p.labels, p.conn);

end % end of regionadjacency