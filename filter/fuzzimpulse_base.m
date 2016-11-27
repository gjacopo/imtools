%% FUZZIMPULSE_BASE - Base function for FUZZIMPULSE.
%
%% Syntax
%     [F,isNoise] = FUZZIMPULSE_BASE(I, K, L);
% 
%% See also
% Related:
% <FUZZIMPULSE.html |FUZZIMPULSE|>. 

%% Function implementation
function [F,varargout] = fuzzimpulse_base(I, K, L)

% image size
[X,Y,C] = size(I);                                                     

% prepare the output
F = I;

%% 
% dealing with multispectral images

if C>1
    if nargout==2,  varargout{1} = false(size(I));  end;
    for c=1:C
        [F(:,:,c), n] = fuzzimpulse_base(I(:,:,c), K, L);
        if nargout==2,  varargout{1}(:,:,c) = n;   end;
    end
    return;
end

%% 
% some useful internal variables

% analyzing windows' parameters
winK = 2*K+1;  % window for impulse noise detection
icenterK = 2*K^2+ 2*K + 1; % sub2ind([winK winK], K+1, K+1);
winL = 2*L + 1; % window for fitering
icenterL =  2*L^2+ 2*L + 1; % sub2ind([winL winL], L+1, L+1);

% indexing matrix
ind = padarray(reshape(1:X*Y,[X, Y]),[1 1],'both','replicate');
% indexing function
fskip = @(A, i, j)  A(2+i:1+i+X,2+j:1+j+Y);

% logical predicates used by fuzzy rules
fAND = @(p1,p2)  p1.*p2;  % min(cat(3,p1,p2),[],3);
fOR = @(p1,p2)  max(cat(3,p1,p2),[],3);
fNOT = @(p1)  (1-p1);

% membership function used by fuzzy rules 
mularge = @(x, a, b) ...              % Fig.(1)
    (x>=b) + ((x-a)./(b-a)).*(x>a&x<b);
% membership function used for filtering
musimilar = @(x, p2, p1, q1, q2) ...  % Fig.(4)
    (x>=p1&x<=q1) + ((x-p2)./(p1-p2)).*(x>p2&x<p1) + ...
    (1-((x-q1)./(q2-q1))).*(x>q1&x<q2);

%% 
% * first detection unit

fun = @(x)  sum(abs(x - repmat(x(icenterK,:),[size(x,1) 1]))) / (winK^2-1);
Obs2 = colfilt(I,[winK winK],'sliding',fun);   % Eq.(4)

Obs1 = colfilt(Obs2,[winK winK],'sliding',@mean);   % Eq.(3)

a = colfilt(Obs1,[winK winK],'sliding',@min);  % Eq.(5)
b = 1.2 * a;   

% fuzzy rule 1
rule1 =  @(p1)  p1;

muimpulse = rule1(mularge(abs(Obs1-Obs2),a,b));

%% 
% * second detection unit

%%
% indices of orientation
%
%          ---------------         -------------
%          | NW | N | NE |         | 1 | 4 | 7 |
%          ---------------         -------------
%          | W  |   | E  |    =>   | 2 | 5 | 8 |
%          ---------------         -------------
%          | SW | S | SE |         | 3 | 6 | 9 |
%          ---------------         -------------
%
% indices of gradients and coordinates skip associated 'related' gradients
NW = struct('ind',1, 'indrel',[+1 -1; -1 +1]);
N =  struct('ind',4, 'indrel',[ 0 -1;  0 +1]);
NE = struct('ind',7, 'indrel',[-1 -1; +1 +1]);
E =  struct('ind',8, 'indrel',[-1  0; +1 0]);
SE = struct('ind',9, 'indrel',[-1 +1; +1 -1]);
S =  struct('ind',6, 'indrel',[ 0 -1;  0 +1]);
SW = struct('ind',3, 'indrel',[-1 -1; +1 +1]);
W =  struct('ind',2, 'indrel',[-1  0; +1 0]);
D0 = 5;
dir_indices = {NW, N, NE, E, SE, S, SW, W, D0};

%%
% compute directional gradients
I = padarray(I,[1 1],'both','replicate');

Delta = cell(9,1);
k = 1;
for j=-1:1  % see above: NW, W, SW, N, central, S, NE, E, SE in this order
    for i=-1:1
        Delta{k} = fskip(I,i,j);      % note: Delta{D0}=I
        k = k+1;
    end
end

for k=1:9
    if k~=5,  Delta{k} = Delta{k} - Delta{D0};   end
end

I = I(2:end-1,2:end-1);

%%
% logical detection variables
gammaimpulse = zeros(X,Y);
gammafree = zeros(X,Y);

a = (Obs1*winK^2 - Obs2) / (winK^2-1);   % Eq.(7)
%a = (colfilt(Obs2,[winK winK],'sliding',@sum) - Obs2) / (winK^2-1); 
b = 1.2 * a;

%%
% composed predicates used in rules 2 and 3
YYY = @(p1,p2,p3)  fAND(fAND(p1,p2),p3);
NNN = @(p1,p2,p3)  fAND(fAND(fNOT(p1),fNOT(p2)),fNOT(p3));
YYN = @(p1,p2,p3)  fAND(fAND(p1,p2),fNOT(p3));
YNN = @(p1,p2,p3)  fAND(fAND(p1,fNOT(p2)),fNOT(p3));

%%
% fuzzy rule 2
rule2 = @(p1,p2,p3)  fOR(YYN(p2,p3,p1), fOR(YNN(p1,p2,p3), fOR(YYN(p1,p2,p3),YYN(p1,p3,p2))));

%%
% fuzzy rule 1
rule3 = @(p1,p2,p3)  fOR(YYY(p1,p2,p3),NNN(p1,p2,p3));

frelgrad = @(D) Delta{D.ind}(fskip(ind,D.indrel(1),D.indrel(2)));
fdirgrad = @(D) Delta{D.ind};

%%
% apply rule
for i=1:8
    % direction 
    D = dir_indices{i};
    % predicates
    p1 = mularge(fdirgrad(D),a,b);
    p2 = mularge(frelgrad(D),a,b);
    p3 = mularge(frelgrad(D),a,b);
    % fuzzy ouput
    gammaimpulse = gammaimpulse + rule2(p1,p2,p3);
    gammafree = gammafree + rule3(p1,p2,p3);
end

isNoise = muimpulse>0 & gammaimpulse>gammafree; % Eq.(8)
nNoise = sum(isNoise(:)==true);

%% 
% filtering

I = padarray(I,[L L], 'both', 'replicate');
% isNoise = padarray(isNoise, [L L], false, 'both');
% xk = find(isNoise);
xk = find(padarray(isNoise, [L L], false, 'both'));

nei = neiposkernel(L, X);
nei(icenterL) = []; % central pixel: inserted later as it corresponds to xk
nei = transpose(nei(:));

xkd = [I(xk) I(repmat(xk,[1 numel(nei)]) + repmat(nei,[nNoise 1]))];

% concatenate and sort the rows (along dimension 2)
xk = sort(xkd, 2);
xij = xk(:,icenterL); % median values in the winL window

rho = sum(diff(xk,1,2),2) / (winL^2-1);    % Eq.(10)

p1 = repmat(xij - rho, [1 winL^2]);
p2 = repmat(xij - 1.1*rho, [1 winL^2]);
q1 = repmat(xij + rho, [1 winL^2]);
q2 = repmat(xij + 1.1*rho, [1 winL^2]);   % Eq.(11)

w = musimilar(xkd, p2, p1, q1, q2);
lambda = w(:,1); % choice of lambda

f = sum(w .* xkd, 2) ./ sum(w, 2);

% F = I(L+1:end-L,L+1:end-L);
% isNoise = isNoise(L+1:end-L,L+1:end-L);
F(isNoise) = lambda .* F(isNoise)  + (1 - lambda) .* f;

if nargout==2,  varargout{1} = isNoise;  end

end % end of fuzzimpulse_base
