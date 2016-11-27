%% ANOMALY_BASE - Perform anomaly detection in multichannel images.
%
%% Description
% 
%
%% Algorithm
% #
%
%% Syntax
%       [] = ();
%
%% Inputs
% *||* : 
%     
%
%% Property [propertyname  propertyvalues]
% *|'|* : 
%     
%
%% Outputs
% *||* : 
%     
%
%% Remarks
% * 
% 
%
%% References
% [TP11]  J. Theiler & L. Prasad: "Overlapping image segmentation for context- 
%      dependent anomaly detection", Proc. SPIE Algorithms and Technologies 
%      for Multispectral, Hyperspectral, and Ultraspectral Imagery 8048, pp.
%      804807-17, 2011.
%      <http://spiedigitallibrary.org/proceedings/resource/2/psisdg/8048/1/804807_1>
%
%% See also
% Related:
% <.html ||>
% <../..//html/.html ||>
% <matlab:webpub(whichpath('')) ||>
% Called:
% <.html ||>
% <../..//html/.html ||>
% <matlab:webpub(whichpath('')) ||>

%% Function implementation
%--------------------------------------------------------------------------
function A = anomaly_base(A, S)
S = compressrange(S); 
A = algorithm4(A, 20);
return
A = algorithm1(A, S);
A = algorithm3(A, S);
end 


%% Subfunctions

%%
% |ALGORITHM1| - Perform anomaly detection over a segmented image; see [TP11],
% p.804808.
%--------------------------------------------------------------------------
function A = algorithm1(A, S)

[X,Y,C] = size(A);
I = 1:X*Y;
CC = accumarrayset(S(:), I(:), false);

A = reshape(A, [X*Y, C]);

%%
% compute mean |M| and covariance |V| of the pixels in the each segment? 
M = cellfun(@(cc) mean(A(cc,:)), CC, 'Uniform', false);
V = cellfun(@(cc) cov(A(cc,:)), CC, 'Uniform', false);

%%
% compute anomalousness of each pixel in the various segment 
M = cellfun(@(cc,m,v) ...
    cellfun(@(a) anomalousness(a, m, v), num2cell(A(cc,:),2)), ...
    CC, M, V, 'Uniform', false);

M = cell2mat(M);
CC = cell2mat(CC);

%%
% output is anomalousness image |A|
A = Inf(X,Y);
A(CC) = M;
end


%%
% |ALGORITHM2| - Double-window anomaly detection; see [TP11], p.804810.
%--------------------------------------------------------------------------
function A = algorithm2(A, S1, S2)

[X,Y,~] = size(A);

%%
% compute |A1| using segmentation |S1| in |ALGORITHM1| 
A1 = algorithm1(A, S1);

%%
% compute |A2| using segmentation |S2| in |ALGORITHM1| 
A = algorithm1(A, S2);

%%
% compute the local pixel min
A = min(A(:),A1(:));
A = reshape(A, [X Y]);

end


%%
% |ALGORITHM3| - Anomaly detection with dilated segments; see [TP11], p.804810.
%--------------------------------------------------------------------------
function Aa = algorithm3(A, S, d, conn)

if isempty(ver('images'))
    warning('algorithm3:missinglibrary', ...
        'Image Processing library not available - use ALGORITHM1 instead');
    Aa = algorithm1(A, S);
    return;
end

if nargin<4,  conn=8;
    if nargin<3,  d = 1;  end
end

[X,Y,C] = size(A);
A = reshape(A, [X*Y, C]);

lS = unique(S); % set of labels

Aa = Inf(X,Y);

if conn==4,      SE = strel('diamond',d);
elseif conn==8,  SE = strel('square',(2*d+1));
end

for l = lS'
    % dilate each labeled region |l| in order to use it as a mask on the input
    % label image: this extracts the original region plus a |d| pixel
    % wide section of any adjacent regions
    R = logical(imdilate(S==l, SE));  R = R(:);
    % compute anomalousness of each pixel in the dilated segments
    m = mean(A(R,:));
    v = cov(A(R,:));
    Aa(R) = min([anomalousness(A(R,:), m, v), Aa(R)], [], 2);
end

end


%%
% |ALGORITHM4| - Anomaly detection with moving window; see [TP11], p.804811.
%--------------------------------------------------------------------------
function A = algorithm4(A, w)
[X,Y,C] = size(A);
A = reshape(A, [X*Y, C]);

%%
% subtract mean from image; this is not strictly necessary, but is numerically
% advantageous
A = A - repmat(mean(A), [X*Y 1]);

%%
% compute |R=a'*a|
Z = cellfun(@(a) mat2rc(transpose(a)*a,'r'), num2cell(A,2), 'Uniform', false);

% %%
% % pack the |C(C+1)/2| distinct components of |R| into |Z|: less memory  
% % but more complexity
% t = triu(true(C,C)); 
% Z = cellfun(@(z) mat2rc(z(t),'r'), Z, 'Uniform', false);

%% 
% transform into a matrix for further filtering
Z = cell2mat(Z); % should be of size (X*Y,C*(C+1)/2)

A = reshape(A, [X Y C]);
Z = reshape(Z, [X Y C*C]);

if ~isempty(ver('images'))
    h = fspecial('average',(2*w+1));
    Z = imfilter(Z,h);
    A = A - imfilter(A,h);
else
    error('algorithm4:missinglibrary', 'not implemented yet');
end

% %%
% % unpack components of |Z| into symmetric matrix |R|

%%
% transform the cells of |Z| back into symmetric covariance matrices
Z = num2cell(reshape(Z, [X*Y, C*C]),2);
Z = cellfun(@(z) reshape(z,[C C]), Z, 'Uniform', false);

A = num2cell(reshape(A, [X*Y, C]), 2);

% the Mahalanobis distance to local mean is anomalousness
Z = cellfun(@(r,a) r - transpose(a)*a, Z, A, 'Uniform', false);
A = cellfun(@(a,r) anomalousness(a, 0, r), A, Z); 
A = reshape(A, [X Y]);
end


%%
% |ANOMALOUSNESS| - Given the mean and the variance of a set of observations,
% compute the anomalousness of those observations. 
%--------------------------------------------------------------------------
function A = anomalousness(a, m, v)
A = (a-m)*(v\transpose(a-m));
end


%%
% |ALGORITHM5| - Simulating one band of a fractal square image of size
% |2^n|.
%--------------------------------------------------------------------------
function A = algorithm5(n, gamma, theta)
A = zeros(2^(n+1), 2^(n+1));

I = 1:n+1;

%%
% width of square; there will be |4^i| of these squares, corresponding to
% the |i|-th level of a quad-tree
W = 2.^(n+1-I);

%%
% define the column (and row) indices for varying |W|
Rj = arrayfun(@(w) ...
    arrayfun(@(i) (1+(i-1)*w:i*w), I, 'Uniform', false), ...
    W, 'Uniform', false); % Rj is a cell of cell
% Rk = Rj;
% (Rj,Rk) is a range that covers a square of width w

%%
% define the ranges for varying |W|
combs = cellfun(@(r) num2cell(allcombs(1:numel(r),1:numel(r)),2), Rj, 'Uniform', false);
range = cellfun(@(r,c) ...
    cellfun(@(i) [mat2rc(r(i(1)),'c') mat2rc(r(i(2)),'c')], combs, 'Uniform', false), ...
    Rj, combs, 'Uniform', false);



%%
% let r be a random scalar
if isempty(ver('images')),  r = rand();
else                        r = normrnd(0,1); % we use a Gaussian distribution
end

%%
% add random value to square of width w


end