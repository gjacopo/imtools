%% IXNEIGHBOURS - Neighbour position indexing.
%
%% Description
% Return the indices of neighbour cells in a matrix. 
%
%% Syntax
%      [ic,icd] = IXNEIGHBOURS(A)
%      [ic,icd] = IXNEIGHBOURS(A, ix);
%      [ic,icd] = IXNEIGHBOURS(A, I);
%      [ic,icd] = IXNEIGHBOURS(A, [], conn);
%
%% Inputs
% *|A|* : |(n,m)| matrix.
%
% *|ix|* : (optional) index vector. 
% 
% *|I|* : (optional) logical matrix of same size as |A|.
%
% *|conn|* : (optional) connectivity index.
%
%% Outputs
% *|ic, icd|* : column vectors of same length where |ic| are the indices of
%     cells in |A| and |icd| are the indices of the neighbour cells.
%
% * |IXNEIGHBOURS(A)| returns all neighbours of all cells in |A|;
% * |IXNEIGHBOURS(A,ix)| returns all neighbours of the cells in index 
%          vector |ix|;
% * |IXNEIGHBOURS(A,I)| returns all neighbours of the cells in the logical
%          matrix |I| that are |true|.
%
%% Remark
% |IXNEIGHBOURS| handles |NaN| values. Hence, it discards cells in |A| that
% are |NaN| both in |ic| and |icd|.
%
%% Examples
%    A = magic(4);
%    A(2,2) = NaN
%    [ic,icd] = ixneighbours(A,3)
%    % construct a sparse adjacency matrix S
%    A = peaks(100);
%    A(A<0) = NaN;
%    nrc = numel(A);
%    [ic,icd] = ixneighbours(A);
%    S = sparse(ic,icd,ones(numel(icd),1),nrc,nrc);
%    spy(S)
% 
%% Acknowledgment
% Original source code by W. Schwanghart available at:
%    http://www.mathworks.com/matlabcentral/fileexchange/14504

%% Function implementation
function [ic,icd] = ixneighbours(X, varargin)

%%
% handle input and error checking
error(nargchk(1, 3, nargin, 'struct'));

siz = size(X);
nrc = siz(1)*siz(2);
In  = isnan(X);

if nargin==1
    X = logical(X);
    if islogical(X)
        method = 'getall';
        nhood  = 8;
    end
    
elseif nargin==2 || nargin==3
    ix = varargin{1};
    if isempty(ix)
        method = 'getall';
    else
        method = 'getsome';
        if islogical(ix)
            if ~isequal(size(X),size(ix))
                error('ixneighbours:inputerror', ...
                    'if I is logical I and X must have same size')
            end
        else
            ixvec = ix(:);
            ix = false(siz);
            ix(ixvec) = true;
        end
        ix = ~In & ix;
    end
        
    if nargin==3
        nhood = varargin{2};
        if ~ismember(nhood(1),[4 8]);
        else
            nhood = nhood(1);
        end
    else
        nhood = 8;
    end
    
end

% replace values in X by index vector
X = reshape((1:nrc)',siz);
X(In) = NaN;    

% pad array
ic  = nan(siz(1)+2,siz(2)+2);
ic(2:end-1,2:end-1) = X;

switch method
    case 'getall'
        I   = ~isnan(ic);
    case 'getsome'
        % Pad logical array
        I = false(siz(1)+2,siz(2)+2);
        I(2:end-1,2:end-1) = ix;
end

icd = zeros(nnz(I),nhood);

% shift logical matrix I across the neighbours
icd(:,1) = ic(I(:,[end 1:end-1]));                % shift to the right                    
icd(:,2) = ic(I([end 1:end-1],:));                % shift down       
icd(:,3) = ic(I(:,[2:end 1]));                    % shift left
icd(:,4) = ic(I([2:end 1],:));                    % shift up

if nhood==8  
    icd(:,5) = ic(I([2:end 1],[end 1:end-1]));        % shift up and right
    icd(:,6) = ic(I([2:end 1],[2:end 1]));            % shift up and left
    icd(:,7) = ic(I([end 1:end-1],[end 1:end-1]));    % shift down and right
    icd(:,8) = ic(I([end 1:end-1],[2:end 1]));        % shift down and left
end

% create output
ic = repmat(ic(I(:)),nhood,1);
icd = icd(:);

% remove NaNs in neighbours
i = isnan(icd);
ic(i) = [];
icd(i) = [];
end % end of ixneighbours              
        
