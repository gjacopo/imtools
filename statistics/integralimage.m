%% INTEGRALIMAGE - Integral image of an image.
%
%% Description
% Compute the integral image of an image defined in [VJ01].
%
% This function computes an integral image so that the value of |F(r,c)| 
% equals |sum(sum(I(1:r,1:c))|. 
% This is nothing else than: |CUMSUM(CUMSUM(I,2),1)|.
%
%% Syntax
%     F =  INTEGRALIMAGE(I);
%
%% References
% [Crow84]  F.C. Crow: "Summed-area tables for texture mapping", SIGGRAPH
%      pp. 207-212, 1984.
%      <http://www.soe.ucsc.edu/classes/cmps160/Fall05/papers/p207-crow.pdf>
%
% [VJ01]  P. Viola and M. Jones, "Robust real-time face detection", 
%      International  Journal of Computer Vision, 57(2):137-154, 2004.
%      <http://research.microsoft.com/~viola/Pubs/Detect/violaJones_IJCV.pdf>
%
%% See also
% Related:
% <INTEGRALHISTO1D.html |INTEGRALHISTO1D|>,
% <INTEGRALHISTO2D.html |INTEGRALHISTO2D|>,
% <INTEGRALHISTOJOINT.html |INTEGRALHISTOJOINT|>.
% Called: 
% <matlab:webpub(whichpath('CUMSUM')) |CUMSUM|>.

%% Function implementation
function F =  integralimage(I)

F = cumsum(cumsum(I,2),1);

end % end of integralimage