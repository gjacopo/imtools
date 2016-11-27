function U = pdelaplace(I,mask)

[M N] = size(I);

% generate and plot the vectorizing indices
G = reshape(1:M*N, [M N]);
G = padarray(G, [1 1], 0,'both');

% generate and display the discrete Laplacian
lapl =  buildlaplacian(G); % sparse(delsq(G));

% modifying Laplacian operator matrix for boundary Neumann conditions
b = reshape((I.*mask)', [M*N 1]);
A = boundlaplacian(lapl, b, M, N);

% solve the linear system
U = reshape(A\b, N, M)';
% [L, U] = luinc(A, 1e-2);
% U = bicgstab(A, b, 1e-3, 100, L, U, b);

end

%--------------------------------------------------------------------------
% different ways of building the Laplacian matrix
function lapl = buildlaplacian(m, varargin)

G = []; % default
if nargin==1 
    if nb_dims(m)~=1,     G = m;
    else                  n = m;   % dimensions passing a square matrix
    end
elseif nargin==2,
    n = varargin{1};   % dimensions passing a rectangular matrix 
end

if ~isempty(G)
    lapl = sparse(delsq(G));

else
    
    % create a partial Laplacian for one "row" of the grid approximation
    % doesn't work for grid widths less than two)
    z = sparse(zeros(1,n));
    z(1) = 4;    z(2) = -1;
    tridiag = kron(speye(m), toeplitz(z));
    
    % now generate the cell prototype for the offdiagonal elements
    z = sparse(zeros(1,m));    z(2) = -1;
    offdiag = toeplitz(z);
    
    %create the discrete Laplacian
    lapl = tridiag + kron(offdiag, speye(n));
end
end


%--------------------------------------------------------------------------
function [A] = boundlaplacian(A, b, varargin)

% If matrix dimensions are passed in, then use them to force boundary
% conditions on the outer rows and columns (ie. the surrounding rectangle)  
% Otherwise, boundary conditions are assumed only for non-zero entries
% If passed in args do not match dimensions for b, they are ignored.

sz = length(b);
m = floor(sqrt(sz));
n = m;
if nargin > 3
  % Rectangular matrix
  arg1 = varargin{1};
  arg2 = varargin{2};
  rows = arg1(1);
  cols = arg2(1);
  if ((rows*cols) == sz)
    m = rows;
    n = cols;
  end

  % now reshape b into a matrix, and force the boundary columns to
  % non-zero
  c = reshape(b, n ,m)';
  c(1, 1:n) = 1;
  c(m, 1:n) = 1;
  c(1:m, 1) = 1;
  c(1:m, n) = 1;
  b = reshape(c', m*n, 1);
end

% Create an identity matrix with ones in the boundary rows, and
% zeros elsewhere
i = find(b);
I = speye(sz);
if length(i) < 1
  A = I;
else
  % This takes FOREVER on large matrices, even sparse ones!
  %A(i,:) = I(i,:);

  % Use a different method; extract diagonals and operate on them only
  % Extend A on left and right so index errors don't occur
  A = [sparse(sz,n) A sparse(sz,n)];

  [B, d] = spdiags(A);
  B(i,:) = 0;
  B(i,3) = 1;
  A = spdiags(B, d, sz, (sz+2*n));

  % Reduce A back to original size
  A = A(:, ((n+1):(sz + n)));
end
end
