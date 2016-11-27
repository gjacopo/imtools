function [SH,TH,pl, pch] = ...
    fuzzrampimpulse_base(I, niter, pilot, method, grd, gthres, kern, ...
    fac, fE, pctch, trans)

error(nargchk(1, 11, nargin, 'struct'));
error(nargoutchk(1, 4, nargout, 'struct'));

if nargin<11,  trans = 'indice';
    if nargin<10,  pctch = 0.01;
        %        if nargin<12,  crit = 'one';   % check only if assertion is true on single band
        if nargin<9,  fE     = 1;       % no correction factor
            if nargin<8,  fac = 'f0';    % factor used in shifting
                if nargin<7,  kern = 'leu';    % original weights for computing the indices
                    if nargin<6,  gthres    = 0;
                        if nargin<5,  grd    = 'sob';   % sobel gradient for estimating the derivatives
                            if nargin<4,  method = 'shift'; % just shift the value, without controlling its extent
                                if nargin<3,  pilot = 'ave';
                                    if nargin<2,   niter = 1;   end
                                end
                            end
                        end
                    end
                end
            end
        end
        %    end
    end
end


%% Parsing parameters
if ~isnumeric(I)
    error('fuzzrampimpulse_base:errorinput','a matrix is required in input');
end


% check the dimension of the input image:
nbdims = nb_dims(I); % instead of ndims
if nbdims<2 || nbdims>3
    error('fuzzrampimpulse_base:inputerror', ...
        'matrix or array of matrices are expected as inputs');
end

% number of spectral components
C = size(I,3);

% possibly overwrite

if C~=3 && strcmp(pilot,'int')                                       
    warning('fuzzrampimpulse_base:inputparameter',...
        'pilot image as intensity implemented only for RGB images');
   pilot = 'ave';
end

%% internal parameters (do not modify)
% Initializing variables

pchange = pctch +1;

% We adopt the following index representation for the directed components
% used for the intensity and gradient indices used (matlab indexing)
%          ---------------         -------------
%          | NW | N | NE |         | 1 | 4 | 7 |
%          ---------------         -------------
%          | W  |   | E  |    =>   | 2 | 5 | 8 |
%          ---------------         -------------
%          | SW | S | SE |         | 3 | 6 | 9 |
%          ---------------         -------------
direction_indices = struct('east', 8, 'northeast', 7, 'north', 4, ...
    'northwest', 1, 'west', 2, 'southwest', 3, ...
    'south', 6, 'southeast', 9, 'central', 5);                         %#ok
%directions = fieldnames(direction_indices);
%ndirections = length(directions);
% note: not used in the following
level_indices = struct( 'L', 1, 'M', 2, 'H', 3);
L = level_indices.('L'); M = level_indices.('M'); H = level_indices.('H');

% Remember: mapping from indexes to subscripts for a 3x3 matrix in matlab
%          -------------------         -------------
%          | 1,1 | 1,2 | 1,3 |         | 1 | 4 | 7 |
%          -------------------         -------------
%          | 2,1 | 2,2 | 2,3 |    =>   | 2 | 5 | 8 |
%          -------------------         -------------
%          | 3,1 | 3,2 | 3,3 |         | 3 | 6 | 9 |
%          -------------------         -------------
%                (ii,jj)          =>    ii+3*(jj-1)

%% construct beforehand the predefined local 3x3 kernels defined for each 
% different zone and level, and used for the estimation of the gradient and
% intensity indices (used once at the beginning of the code).
% Intensity and gradient masks are built for both zones 1 and 2, and all
% the other zones derived by rotation
switch kern
    case 'leu'
        matI = local3x3kernel('ker','i0','norm',true);
        matG = local3x3kernel('ker','g0','norm',true);
    case 'new'
        matI = local3x3kernel('ker','i1','norm',true);
        matG = local3x3kernel('ker','g1','norm',true);
end
% mote: matI and matG are indexed by [size(x,y),zone,level]


%% Main computation through iterative filtering

% % possibly resize the input matrix
% if interp
%     A = upscalexy(I,[2 2],'cubic');
% else
     A = I;
% end
% initialize the output matrix
% SH = A;

% dimension of the frame
[X,Y] = size(A(:,:,1)); 
XY = X * Y; % numel(A(:,:,1));

% index of all pixels in the input image
% pixindex = reshape(1:XY,[X Y]);
% indexes of the border pixels
% pixbord = [1:X, (1:Y-2)*X+1, (2:Y-1)*X, (1+X*(Y-1)):XY]'; 

% create the 'pilot' for gradient orientation
if C==3 && strcmp(pilot,'bright')
   pilot = rgb2gray(A); % pilot will be the brightness image
elseif C>=2
    pilot = sum(A,3) / C; % pilot will be the average image
end

% construct the variation sparse matrices measuring the amount of change
% occurring in the image after each of the iterative filtering
deltaI = zeros(XY,C);   
dirdeltaI = deltaI; 

% index of transition pixels for each step of the iterations
Itrans = cell(niter);

% set the number of estimated gradient
nC = C + (C>1); % ie: ng=1 if Z==1, ng=Z+1 otherwise

G = zeros(X,Y,nC);

pl = A(110, :);
pch = [];

TH = zeros(X,Y);

% proceed iteratively
for iter=1:niter
        
    
    %% Computation of the directional derivatives through Gaussian smoothing and
    % differentiation
    % compute the gradient (gy: vertical, gx: horizontal) for each channel
    [gy,gx] = grdmask_base(A, grd, 'ij');
    % same as computing first: [gx,gy]=grdsmooth(I,sigma,p.der,hsize,'xy');
    % and then take the vector orthogonal to the gradient: tmp=gx; gx=gy; gy=-tmp;
    gy = -gy;   
    
    % note that the output directional derivatives have size [X Y C]
    if C>1
        % norm channel by channel
        G(:,:,1:C) = sqrt(gx.^2 + gy.^2);
        % update the value of the gradient to set it to the gradient of
        % the average image
        gx = sum(gx,3) / C;
        gy = sum(gy,3) / C;
    end
   % Theta = mod(atan2(gy,gx),pi);
   Theta = atan2(gy,gx);
    figure, imagesc(Theta), colorbar
    G(:,:,nC) = mean(G(:,:,1:C),3);
    figure, imagesc(rescale(G(:,:,nC))), colormap gray

     % find the orientation and the interpolation parameters over the image
    [Zones,Omega] = localorientzone(Theta,8);
    % compute the compensation factor
    S = 1 - (1-sqrt(2.)) * Omega;

    
    %% local estimation of intensity indices
    
    % prior computation of the intensity indices over the different
    % spectral components
    mI = localorientfeature(A, 'filt', 'mean', ... % 'filt','med'
        'Kernel',matI,'Zones',Zones,'Omega',Omega);   
    % mI = round(mI);

   
    %% local estimation and characterization of ramp/transition pixels
    
    Iramp = maptransition(mI,trans,'const','strong');
   figure, imagesc(reshape(Iramp,[X,Y])), colormap gray;
   
    % get rid of flat area:
    if gthres > 0
        m = max(max(G(:,:,nC)));
    else 
        m=0;
    end
    Iramp = Iramp & reshape(G(:,:,nC),[XY 1])>m*gthres;
            
    % if no consideration for this condition:  Iramp = ones(XY,1);
    %figure, imagesc(reshape(Iramp(:,:,1),X,Y)), axis image, colormap gray;
    Iramp = find(Iramp);
    
    % proceed only if such pixels have been found
    if isempty(Iramp) % this is very improbable
        break;
    end
    
    %% local estimation of the gradient indices
    %mG = zeros(lenght(Iramp), nelevels, Z);
    mG  = localorientfeature(G,'filt','mean', ...
        'Kernel',matG,'Zones',Zones,'Omega',Omega);
    
    % reduce the problem to potential ramp pixels: restrict the set of
    % pixels which are examined to pixels on the ramp
    mI = mI(Iramp,:,:);
    mG = mG(Iramp,:,:);
    
    if strcmp(method,'control')
        D = deltaI(Iramp,:);
        ID = dirdeltaI(Iramp,:);
    else
        D = zeros(length(Iramp),C);
        ID = [];
    end
    S = S(Iramp);
    
    % extraction of transition pixels (located on a ramp) with the
    % criterion GM>GH and GM>GL
    iR = mG(:,M,nC)>mG(:,H,nC) & mG(:,M,nC)>mG(:,L,nC);

    Itrans{iter} = Iramp(iR); % a subset of the ramp pixels
    a = zeros(X,Y);  a(Itrans{iter}) = 1;
   figure, imagesc(a),colormap gray, title('ramp')
    
   
    %% Image sharpening
    
    % reshape the input and initialize the output
    A = reshape(A,[XY C]);
    SH = A;
    
    % extract ramp pixels of type 1:    GL <(or<=) GM <= GH
    iR = findcase(method, mG(:,H,nC), mG(:,L,nC), mG(:,M,nC));
    % note : iR are the  coordinates of the pixels considered in the domain
    % of the reduced image and Iramp(iR) are the correponding coordinates in
    % the domain of the original image
    % possibly update those pixels
    if ~isempty(iR)
        for c=1:C
            F = factorvalue(fac, mG(iR,H,c), mG(iR,L,c), mG(iR,M,c));
           if strcmp(method,'shift')
                R = adjustleu(F, mI(iR,L,c), mI(iR,M,c), S(iR), fE);
                SH(Iramp(iR),c) = updateleu(A(Iramp(iR),c), R, -1);
            elseif strcmp(method,'control')
                R = adjustcontrol(F, mI(iR,L,c), mI(iR,M,c), D(iR,c), ID(iR,c), ...
                    S(iR), fE);
                 SH(Iramp(iR),c) = updatecontrol(A(Iramp(iR),c),R,mI(iR,L,c));
            end
        end
    end
    
    % extract ramp pixels of type 2:    GH <(or<=) GM <= GL
    iR = findcase(method, mG(:,L,nC), mG(:,H,nC), mG(:,M,nC));
    if ~isempty(iR)
        for c=1:C
            F = factorvalue(fac, mG(iR,L,c), mG(iR,H,c), mG(iR,M,c));
            if strcmp(method,'shift')
                R = adjustleu(F, mI(iR,H,c), mI(iR,M,c), S(iR), fE);
                SH(Iramp(iR),c) = updateleu(A(Iramp(iR),c), R, 1);
            elseif strcmp(method,'control')
                R = adjustcontrol(F, mI(iR,H,c), mI(iR,M,c), D(iR,c), ID(iR,c), ...
                    S(iR), fE);                 
                SH(Iramp(iR),c) = updatecontrol(A(Iramp(iR),c),R,mI(iR,H,c));
            end
        end
    end
    
    % SH = round(SH);
    
    if niter>1
       % updates:
        %  - matrices delta of intensity variation changes
        delta = A(Iramp,:) - SH(Iramp,:);        
        %  - matrices dirdelta of change in intensity variation direction
        for c=1:C
            dirdeltaI(Iramp(delta(:,c) .* deltaI(Iramp,c) < 0),c) = 1; % sign change
        end
        deltaI(Iramp,:) = delta;
        pchange = sum(abs(delta),2)>eps; % matrix of modified pixels
        pchange = sum(pchange(:)) / XY;  % pct of change
        % note : the first sum: sum(abs(delta),3) operates over the channnel
        % account for pixels modified in any of their channel
        
    end
        
    %% Process for update for next loop in the iteration
    A = reshape(SH, [X, Y, C]);
    
     pl = [pl ; A(110, :)];
     pch = [pch ; pchange];
    
    if pchange <= pctch
        pchange
        pctch
        break;
    end

% TH: final ramp after the last iteration
iR = findramp(method, mG(:,L,nC), mG(:,H,nC), mG(:,M,nC));
TH(Iramp(iR)) = TH(Iramp(iR))+1;
    
end


% final output
SH = A;

end% end of rampsharp


%% adjustment estimation

% -------------------------------------------------------------------------
function iR = findcase(method, G1, G2, Gm)
iR = G1>=Gm & Gm>G2; % standard Leu condition
if ~strcmp(method,'leu')
    iR = iR | (G1>=Gm & Gm==G2); % add flexible condition
end
%iR = find(iR);
end
% end of findcase

% -------------------------------------------------------------------------
function iR = findramp(method, G1, G2, Gm)
iR = Gm>=G1 & Gm>=G2; % 
%iR = find(iR);
end
% end of findramp

% -------------------------------------------------------------------------
function SH = updatecontrol(A, R, I1)
SH =  (A > I1) .* max(A - R, I1) + (A <= I1) .* min(A + R, I1);
end
% end of updatecase


% -------------------------------------------------------------------------
function SH = updateleu(A, R, s)
SH =  A + s * R;
end
% end of updateleu


% -------------------------------------------------------------------------
function R = adjustcontrol(F, I1, I2, D, ID, S, fE)
R = 0.5 + fE .* F .* S .* abs(I1-I2);
iR0 = D~=0 & ID==1;
% control the current correction amount by the previous correction amount
if ~isempty(iR0)
    R(iR0) = max(min(R(iR0), abs(D(iR0))-1), 0);
end
end
% end of adjustcase


% -------------------------------------------------------------------------
function R = adjustleu(F, I1, I2, S, fE)
R = fE .* S .* abs(I1-I2);
R = (F >= 0.5) .* R + 2. * (F < 0.5) .* F .* R;
end
% end of adjustleu


% -------------------------------------------------------------------------
function F = factorvalue(alpha, G1, G2, Gm)
% compute the correction factor based on the estimated gradient indices
% Gm, Gl (either G1 or G2) and Gh (ibid)
% G1 and G2 stand for Gl or Gh depending on the part of the ramp (lower or
% higher) the pixel belongs to
switch alpha
    case {'leu','f0'}  % original leu
        F = (G1 - Gm) ./ (Gm - G2); % leu-like
    case 'f1'  % default choice
        F = (G2 + Gm) ./ (1 + 2*Gm);
    case 'f2'
        F = (G1 + G2 - Gm) ./ (1 + 3*Gm);
    case 'f3'
        F = (G1 + G2) ./ (1 + 2*Gm);
    case 'f4'
        F = (G1 + Gm) ./ (1 + 2*Gm);
    case 'f5'   
        F = (G1 + G2 + Gm) ./ (1 + 3*Gm);
end
end
% end of factorvalue


