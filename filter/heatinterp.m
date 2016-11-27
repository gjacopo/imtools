function Mi = heatinterp(I,mask,iter,method,N,K,dt)
% I    source image (2D gray-level matrix) for diffusion
% iter    number of iterations
% dt   time increment (0 < dt <= 0.25, default 0.2)
% method =  'lin':  Linear diffusion (constant c=1).
%           'pm1': perona-malik, c=exp{-(|grad(J)|/K)^2} [PM90]
%           'pm2': perona-malik, c=1/{1+(|grad(J)|/K)^2} [PM90]
%           'rmp': complex valued - ramp preserving [GSZ01]

if ~exist('method','var')
   method='lin';
end
if ~exist('N','var')
   N=1;
end
if ~exist('K','var')
   K=1;
end
if ~exist('dt','var')
   dt=0.2;
end

Mi = I;
Mi(mask==1) = mean(I(mask==0));
for i=1:iter
    Mi(mask==0) = I(mask==0);
    Mi = diffusion(Mi,method,N,K,dt);
end

end

function U = heatneumann(I,mask,iter,method,N,K,dt)

%  grid and time stepping parameters
A = padarray(I, [1 1], 'replicate','both');
[m,n] = size(I);
[m2,n2] = size(A);  
L =   1;  
h = 1 ./ [m n];  
T = .5;  dt = T/N;

x = -h(1)/2:h(1):L+h(1)/2; y = -h(2)/2:h(2):L+h(2)/2;
[xg yg] = meshgrid(x,y); 
tg = 0:dt:T;

%  Diffusion coefficient
D = 1;
lambda = D*dt/(h(1)*h(2));

% generate & plot the vectorizing indices
G = reshape(1:m*n, [m n]);
G = padarray(G, [1 1], 0,'both');
%   spy(G)
%   disp('pause'), disp(' '), pause

% boundaries
dmask = imdilate(mask,strel('square',3))-mask;
imagesc(dmask), colormap gray

% Generate and display the discrete Laplacian, modified for 
% Neumann conditions
lapl = sparse(-delsq(G)/(h(1)*h(2)));
lapl(dmask>0) = lapl(dmask>0) + 1/(h(1)*h(2));


% Number of interior points (count when G_index > 0)   
nintern = sum(G(:)>0);

%  Initial conditions
u = I;
u(mask==0) = I(mask==0);
u(mask==1) = mean(I(mask==0));

U = G;   U(G>0) = full(u(G(G>0)));
%   mesh(U);

rhs = u(:) ;

%  PDE solve with operation count
% Crank-Nicolson, 2d Heat Equation

hmat = speye(nintern) - D * dt * lapl;
% [LL,UU] = luinc(hmat,1e-4);
RR = cholinc(hmat,1e-6);

for k=1:N
    time = k*dt;
    
    % iteration
    %   z = hmat\rhs;
    %   z = pcg(hmat,rhs,1d-14,100);
    %   z = gmres(hmat,rhs,50);

    %   u = cgs(hmat,rhs,1d-8,10,LL,UU);
    u = pcg(hmat,rhs,1d-10,20,RR',RR);
    %  plot the solution (with BCs & corners = 0)
    
    U = G;
    U(G>0) = full(u(G(G>0)));
    
    % update right-hand side
    %  Initial conditions and boundary conditions
    rhs = u ;
    rhs(mask==0) = I(mask==0);
        
end
end