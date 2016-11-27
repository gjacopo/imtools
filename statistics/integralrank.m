%% INTEGRALRANK - Compute rank image.
%
%% Description
% Use an integral histogram based approach to compute local rank filtered
% image.
%
%% Syntax
%     R = INTEGRALRANK(I, rank, win, res);
%
%% See also
% Related:
% <INTEGRALIMAGE.html |INTEGRALIMAGE|>,
% <SLIDEHISTOFUN_BASE.html |SLIDEHISTOFUN_BASE|>.
% Called:
% <INTEGRALHISTO1D.html |INTEGRALHISTO1D|>,
% <matlab:webpub(whichpath('MESHGRID')) |MESHGRID|>,
% <matlab:webpub(whichpath('CUMSUM')) |CUMSUM|>,
% <matlab:webpub(whichpath('UNIQUE')) |UNIQUE|>,
% <matlab:webpub(whichpath('FIND')) |FIND|>.

%% Function implementation
function R = integralrank(I, rank, win, res)

%%
% checking parameters and setting internal variables

[X,Y,C] = size(I);                                                     %#ok

% ensure that the window size is odd
if mod(win,2)==0,  win = win + 1;  end
ws = floor(win / 2); % half size window

%% 
% computation

% pad the input image
A = padarray(I, [ws ws],'replicate','both'); 

IH = integralhisto1d(A, res);
% IH = IH(1+pad:M+pad, 1+pad:N+pad, :);
IH = padarray( IH, [1 1 0], 0, 'both');

M = X+2*(1+ws);
N = Y+2*(1+ws);

IH = reshape(IH,[M*N res]);
size(IH)

[indY,indX] = meshgrid(1:X, 1:Y);
indY = indY(:); indX = indX(:);
I0 = (indX+1+2*ws) + M*(indY+1+2*ws);
max(I0(:))

I1 = indX + M*(indY+1+2*ws);
I2 = (indX+1+2*ws) + M*indY;
I00 = indX + M*indY;

pI = IH(I0,:) - IH(I1,:) - IH(I2,:) + IH(I00,:);
pI = cumsum(pI,2);
[xy,r] = find(pI>=rank);
[~,ind] = unique(xy);

R = reshape(r(ind),[X,Y]);

% R = zeros(X,Y);
% for x=1:X
%     for y=1:Y
%         pI = IH(x+1+2*ws,y+1+2*ws,:) - ...
%             IH(x,y+1+2*ws,:) - IH(x+1+2*ws,y,:) + IH(x,y,:);
%         pI = cumsum(squeeze(pI));
%        % [x y]
%        % find(pI>=rank,1,'first')
%         R(x,y) = find(pI>=rank,1,'first');
%     end
% end

end % end of integralrank