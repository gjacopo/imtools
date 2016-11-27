%% LOCALSUM - Compute local sum of an image in a square window.
%
%% Description
% This function is nothing else than a copy/paste of the subfunction with
% same name in the Matlab built-in function |NORMXCORR2| (implementation of
% Normalized two-dimensional cross-correlation).
% 
%% Syntax
%     S = LOCALSUM(I, m);
%     S = LOCALSUM(I, m, n);
%
%% Inputs
% *|I|* : input image (matrix).
%
% *|m, n|* : dimensions |(m,n)| of the neighbourhood window over which the
%     entries of the input matrix are summed; when only m is passed, a square
%     window |(m,m)| is considered.
% 
%% Output
% *|S|* : local summed matrix.
%
%% Acknowledgment
% from |NORMXCORR2|: The algorithm depends on precomputing running sums as
% described in "Fast Normalized Cross-Correlation", by J. P. Lewis, Industrial
% Light & Magic. 
% See <http://www.idiom.com/~zilla/Papers/nvisionInterface/nip.html>
%
%% See also
% Related:
% <matlab:webpub(whichpath('NORMXCORR2')) |NORMXCORR2|>,
% <matlab:webpub(whichpath('SUM')) |SUM|>.
% Called:
% <matlab:webpub(whichpath('CUMSUM')) |CUMSUM|>.

%% Function implementation
function local_sum_A = localsum(A, m, n)

if nargin==2, n=m;  end

B = padarray(A,[m n]);
s = cumsum(B,1);
c = s(1+m:end-1,:)-s(1:end-m-1,:);
s = cumsum(c,2);
local_sum_A = s(:,1+n:end-1)-s(:,1:end-n-1);

end % end of localsum