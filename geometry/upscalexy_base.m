%% UPSCALEXY_BASE - Base function for UPSCALEXY.
%
%% Syntax
%      upI = UPSCALEXY_BASE(I, S, method);
% 
%% See also 
% Related:
% <UPSCALEXY.html |UPSCALEXY|>, 
% <DOWNSCALEXY_BASE.html |DOWNSCALEXY_BASE|>.
% Called:
% <matlab:webpub(whichpath('UPSAMPLE')) |UPSAMPLE|>,
% <matlab:webpub(whichpath('CUMSUM')) |CUMSUM|>,
% <matlab:webpub(whichpath('MESHGRID')) |MESHGRID|>,
% <matlab:webpub(whichpath('INTERP2')) |INTERP2|>.

%% Function implementation
function upI = upscalexy_base(I, S, method)

% dealing with multichannel images
% if C>1
%     upI = upscalexy_base(I(:,:,1), S, method); 
%     for c=2:C
%         upI = cat(3, upI, upscalexy_base(I(:,:,c), S, method));
%     end
%     return;
% end
% the function UPSCALEXY_BASE runs already in its current version for multichannel
% images

%%
% internal variables

if length(S)==1,        S = [S S 1];                           
elseif length(S)==2,    S = [S 1];
end
% we must have S=[sx sy 1] with size(S)=3

%  size (and number of dimensions) of the input matrix
sI = size(I);
if numel(sI) == 2, sI(3) = 1; end;

if all(S==1) % just don't bother...
    upI = I;    return;
end

%% 
% upscaling

switch method
    
    case 'replicate'
        for i = length(sI):-1:1
            % one index vector into A for each dim
            H = zeros(sI(i)*S(i),1);
            % put ones in correct places
            H(1:S(i):sI(i)*S(i)) = 1;  
            % cumsumming creates the correct order
            T{i} = cumsum(H);  
        end
        % feed the indices into A
        upI = I(T{:}); 
        
    case 'upsample' % upsamples signal by inserting 0 between input samples
        upI = zeros(sI(1)*S(1), sI(2)*S(2),sI(3));
        for c=1:sI(3)
            upI(:,:,c) = upsample(upsample(I(:,:,c),S(1))',S(2))';
        end
        
    otherwise
        % create the output image
        upI = zeros(sI(1)*S(1), sI(2)*S(2),sI(3));
        % create the grid for interpolation
        [x,y] = meshgrid(1:sI(2),1:sI(1));
        [xi,yi] = meshgrid( 1:1/S(2):(sI(2)+1-1/S(2)), ...
                            1:1/S(1):(sI(1)+1-1/S(1)));
        xi(:,end) = sI(2); yi(end,:) = sI(1);
        % 2D interpolation
        for c=1:sI(3)
            upI(:,:,c) = interp2(x, y, I(:,:,c), xi, yi, method);
        end
        
end
end % end of upscalexy_base


