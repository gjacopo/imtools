function O = slidecovfun_base(I, func, win)
%% References
%   [VJ01]  P. Viola and M. Jones, "Robust Real-Time Face Detection", 
%      International  Journal of Computer Vision, 57(2):137-154, 2004.
%      <http://research.microsoft.com/~viola/Pubs/Detect/violaJones_IJCV.pdf>
%
%   [Pori05]  F. Porikli: "Integral histogram: a fast way to extract histogram
%      features", Proc. IEEE CVPR, pp. 829?836, 2005.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1467353>
%
%   [PT06]  F. Porikli and O. Tuzel: "Fast construction of covariance matrices
%      for arbitrary size image windows", Proc. IEEE ICIP, 2006.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4106846>
%
%   [TPM06]  O. Tuzel, F. Porikli and P. Meer: "Region covariance: a fast 
%      descriptor for detection and classification", Proc. ECCV, LNCS 3952, 
%      pp. 589-600, 2006.
%      <http://www.springerlink.com/content/r7007734n986208t/>
% 
%% See also
% Called:
% <INTEGRALIMAGE.html |INTEGRALIMAGE|>,
% <../../vista/algebra/ACLLCOMBS.html |ACLLCOMBS|>,
% <../../vista/algebra/TRIUNPACK.html |TRIUNPACK|>,
% <matlab:webpub(whichpath('NUM2CELL')) |NUM2CELL|>,

%%
% * checking/setting the parameters/data

[X,Y,C] = size(I);
% if nargin<5,  bin = max(I(:)) + 1;   end

% ensure that the window size is odd
if mod(win,2)==0,  win = win + 1;  end
ws = floor(win / 2); % half size window

if strcmpi(class(func),'function_handle'),   func = {func};  end
nf = numel(func);

% prepare the output data
O = cell(nf,1);
for i=1:nf,   O{i} = zeros(X,Y);  end

% * pad the input image
A = padarray(I, [ws+1 ws+1],'symmetric','both');
[M,N] = size(A(:,:,1));

if strcmpi(class(func),'function_handle'),   func = {func};  end
nf = numel(func);

% prepare the output data
O = cell(nf,1);
for i=1:nf,   O{i} = zeros(X,Y);  end

%%
% * compute the matrix (or tensor when |C>1|) of |C| integral image(s), cf.
% Eq.(10) of [TPM06] or Eq.(6) of [PT06]; we transform the integral images
% into a cell indexed by |C| of matrices of integral values

P = integralimage(A); % cf. Eq.(8) of [TPM06] or Eq.(4) of [PT06]

%%
% transform the |C| matrices of |(M,N)| integral images into a matrix of
% dimension |(M*N,C)|
P = reshape(P,[M*N C]); % in the case C=1, it is equivalent to P = P(:);

%%
% finally, reconvert |P| into a cell indexed by the position and storing
% the |C| vectorial integral values
P = num2cell(P,2); % each one of the |M*N| cells contains |C| integral values

%%
% * compute the matrix (ibid, tensor) of the second order integral image(s);
% cf. Eq.(10) of [TPM06] or Eq.(6) of [PT06]
if C>1
    %%
    % compute the order 2 integral images for all (non redundant) combinations
    % of vectorial channels
    Q = cellfun(@(c) integralimage(A(:,:,c(1)) .* A(:,:,c(2))), ...
        num2cell(unique(sort(allcombs(1:C,1:C), 2), 'rows'),2), ...
        'Uniform', false); % cf. Eq.(9) of [TPM06] or Eq.(5) of [PT06]
    % note that ALLCOMBS(1:C,1:C)=FLIPLR(ALLVCOMBS(1:C,2))
    nQ = numel(Q); % in fact, |nQ=C*(C+1)/2|
    
    %%
    % convert the tensor into a a matrix of size |(M*N,C*(C+1)/2)|
    % containing the non redundant entries of 2-combinations of (1:C)
    Q = cell2mat(cellfun(@(x) x(:), Q, 'Uniform', false));
    Q = reshape(Q,[M*N,nQ]);
    
    %%
    % transform back into |M*N| cells contains the |(C,C)| covariance
    % matrices
    Q = num2cell(Q,2); % each one of the |M*N| cells contains |C(C+1)/2| integral values
    Q = cellfun(@(q) triunpack(q,C,'s'), Q, 'Uniform', false);
    
else
    %%
    % simpler |C=1| case
    Q = integralimage(A.^2);
    Q = num2cell(Q(:));
end

% slide
for x=1:X
    for y=1:Y
        wC = covint(P, Q, X, x, y, ws);
        for i=1:nf
            O{i}(x,y) = func{i}(wC);
        end
    end
end

%   %----------------------------------------------------------------------
    function [wh, H] = slidecovint(C, P, Q, ws, x);             %#ok
        % sliding function
    end
%   %----------------------------------------------------------------------
    function wc = covint(P, Q, X, x, y, ws)          
        xp = x + (2*ws+1);  yp = y + (2*ws+1);
        xpyp = (yp-1)*X+xp;  xy = (y-1)*X+x;
        xpy = (y-1)*X+xp;  xyp = (yp-1)*X+x;
        n = (xp-x) * (yp-y);
        p = P{xpyp} + P{xy} - P{xyp} - P{xpy};
        q = Q{xy} + Q{xpyp} - Q{xpy} - Q{xyp};
        wc = 1/(n-1) * (q - 1/n * transpose(p) * p);
    end
%   %----------------------------------------------------------------------

if nf==1,  O = O{1};  end

end % end of slidecovfun_base
