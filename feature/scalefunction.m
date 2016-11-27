%% SCALEFUNCTION - Representation of typical edges and steps.
%
%% Description
% Creates a function 1D or 2D with some typical edges and steps.
%
%% Syntax
%     F = SCALEFUNCTION( sbin, dim, bound );
%     
%% Inputs
% *|sbin|* : the image will have size |11*sbin|; must be a multiple of 3 (or
%     will be forced.
%
% *|dim|* : result will be a graph (|dim=1|) or an image (|dim=2|); default:
%     |dim=1|.
%
% *|bound|* : vector |(3,1)| of the |[min,med,max]| for the values of the 
%     function; default: |bound=[0.8,1,2]|.
%
%% Output
% *|F|* : 1D or 2D-function representing some typical edges and steps.

%% Function implementation
function F = scalefunction( sbin, dim, bound )

if exist('dim','var')~=1,  dim=1; end;
if exist('bound','var')~=1
  Min=0.8; Med = 1.; Max = 1.2;
else 
  Min=bound(1); Med=bound(2); Max=bound(3);
end;

Max2 = Med+ 2.*(Max-Med)/3.;
Min2 = Med - 2.*(Med-Min)/3.;

q=floor(sbin/3.);
if(sbin~=3*q)
  warning('scalefunction:inputwarning', ...
      ['Parameter sbin=', num2str(sbin), 'changed to sbin=' num2str(3*q)] );
end
sbin=3*q;

size = 11 * sbin; 
ncol = size;
% nraw=size ou 1
nraw = (dim==2)*size + (dim==1);

x = 1:size; y = 1:nraw;
[X,~] = meshgrid(x,y);

F = ones(nraw, ncol);

%%%%  Edge type 1
% branche 1
alpha = (Max-Med) / sbin;
beta = Med - alpha * sbin;
F(:, sbin+1:2*sbin) = alpha * X(:, sbin+1:2*sbin) + ...
    beta;
% branche 2
alphap = (Med-Max) / sbin;
betap = Med + 3*alpha * sbin;
F(:, 2*sbin+1:3*sbin) = alphap * X(:, 2*sbin+1:3*sbin) ...
    + betap;

%%%% Gaussian
sigma=2;
fac = exp(-1/(2*sigma));
[x,~] = meshgrid(linspace(-1,1,2*sbin),1:nraw);
F(:,4*sbin+1:6*sbin) = 1 -fac + exp( - x.^2 / (2*sigma));
				   
%%%% Edge type 2
x1=q*23; x2=q*25;
% branche 1
alpha= (Max2-Med) / (2*q);
beta = Med - 7*alpha * sbin;
F(:,7*sbin+1:x1) = alpha * X(:, 7*sbin+1:x1) + beta;
% branche 3
alphap= (Med-Min2) / (2*q);
betap = Med - 9*alphap * sbin;
F(:,x2:9*sbin) = alphap * X(:, x2:9*sbin) + betap;
% branche 2
alpha2 = (Min2-Max2) / (2*q);
beta2 = Max2 - alpha2 * 23 * q;
F(:,x1:x2) = alpha2 * X(:, x1:x2) + beta2;

%%%% Step
F(:,10*sbin+1:ncol) = Min;

if dim==1  
  figure, plot(X,F)
elseif dim==2  
  figure, surf(F); colormap jet;
  figure, imagesc(F), colormap gray, axis image;
end

end
