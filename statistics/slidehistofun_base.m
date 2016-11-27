%% SLIDEHISTOFUN_BASE - Base function for SLIDEHISTOFUN.
%
%% Syntax
%    O = SLIDEHISTOFUN_BASE(I, method, func, win, bin);
%
%% See also
% Related:
% <SLIDEHISTOFUN.html |SLIDEHISTOFUN|>,
% Called:
% <INTEGRALHISTO1D.html |INTEGRALHISTO1D|>,
% <INTEGRALHISTOJOINT.html |INTEGRALHISTOJOINT|>,
% <matlab:webpub(whichpath('HISTC')) |HISTC|>.

%% Function implementation
function O = slidehistofun_base(I, method, func, win, bin)

%%
% checking/setting the parameters/data

[X,Y,C] = size(I);                                                     %#ok

% if nargin<5,  bin = max(I(:)) + 1;   end

% ensure that the window size is odd
if mod(win,2)==0,  win = win + 1;  end
ws = floor(win / 2); % half size window

if strcmpi(class(func),'function_handle'),   func = {func};  end
nf = numel(func);

% prepare the output data
O = cell(nf,1);
for i=1:nf,   O{i} = zeros(X,Y);  end

%%
% initializing function

inithisto = str2func(['inithisto' method]);

%   %----------------------------------------------------------------------
    function H = inithistoint(A, dum1, bin)                            %#ok
        H = integralhisto1d(A, -bin);
        H = padarray( H, [1 1 0], 0, 'both');
    end
%   %----------------------------------------------------------------------
    function H = inithistojoint(A, dum1, bin)                          %#ok
        H = integralhistojoint(A, A, -bin);
        H = padarray( H, [1 1 0], 0, 'both');
    end
%   %----------------------------------------------------------------------
    function cH = inithistodist(A, ws, bin)                            %#ok
        % % initial kernel histogram: count the number of entries accumulated
        % % in each bin
        % H = A(1:2*ws+1,1:2*ws+1);
        % %H = accumarray(H(:), 1, [res 1], @sum);
        % H = hist(H(:),res);
        
        % column histograms
        cH = histc(A(2:2*(ws+1), :), 0:bin-1); % matrix of size (res,size(A,2))
    end
%   %----------------------------------------------------------------------

%%
% sliding function

slidehisto = str2func(['slidehisto' method]);

%   %----------------------------------------------------------------------
    function [wh, H] = slidehistoint(A, H, x, Y, ws, bin)              %#ok             
        wh = zeros(bin,Y);
        for y=1:Y
            wh(:,y) = H(x+1+2*ws,y+1+2*ws,:) - ...
                H(x,y+1+2*ws,:) - H(x+1+2*ws,y,:) + H(x,y,:);
        end
        % do nothing to H: dummy variable here
        % do not use A: dummy variable
    end
%   %----------------------------------------------------------------------
    function [wh, H] = slidehistojoint(A, H, x, Y, ws, bin)            %#ok
        wh = zeros(bin,Y);
        for y=1:Y
            wh(:,y) = H(x+1+2*ws,y+1+2*ws,:) - ...
                H(x,y+1+2*ws,:) - H(x+1+2*ws,y,:) + H(x,y,:);
        end
        % do not use A: dummy variable
    end
%   %----------------------------------------------------------------------
    function [wH, cH] = slidehistodist(A, cH, x, Y, ws, bin)           %#ok
        x = x+1;
        xx = x+2*ws;
        
        if x>1
            % update the (2*ws+1) leftmost column histograms
            i0 = (0:2*ws)*bin+A(x-1,1:1+2*ws)+1;  cH(i0) = cH(i0) - 1;
            i1 = (0:2*ws)*bin+A(xx,1:1+2*ws)+1;   cH(i1) = cH(i1) + 1;
        end
        
        % allocate  the kernel histogram for the current line
        wH = zeros(bin,Y);
        % initialize the kernel histogram for the first pixel/neighbourhood
        % with the recently updated (2*ws+1) leftmost column histograms
        wH(:,1) = sum(cH(:,1:1+2*ws),2);
        
        for y=2:Y
            % update the rightmost column histogram
            yy = y+2*ws;
            i0 = A(x-1,yy)+1;  cH(i0, yy) = cH(i0, yy) - 1;
            i1 = A(xx,yy)+1;   cH(i1, yy) = cH(i1, yy) + 1;
            % update the kernel
            wH(:,y) = wH(:,y-1) + cH(:,yy) - cH(:,y-1);
        end
    end
%   %----------------------------------------------------------------------

%%
% computation

% pad the input image
A = padarray(I, [ws+1 ws+1],'symmetric','both'); 
% we choose 'symmetric' because it makes the calculation of the kernel 
% histogram and the column histograms of the first line of the image much
% easier

H = inithisto(A, ws, bin);

% slide
for x=1:X
    [wH, H] = slidehisto(A, H, x, Y, ws, bin);
    for i=1:nf
        O{i}(x,:) = func{i}(wH);
    end
end

if nf==1,  O = O{1};  end

end % end of slidehistofun_base
