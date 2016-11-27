%% IM2FRONT_BASE - Base function for IM2FRONT.
%
%% Syntax
%    D = IM2FRONT_BASE();
%
%% See also
% Related:
% <IM2FRONT.html |IM2FRONT|>, 
% <FMM_BASE.html |FMM_BASE|>.
% Called: 
% <POTENTIAL2FRONT.html |POTENTIAL2FRONT|>,
% <IM2POTENTIAL.html |IM2POTENTIAL|>.

%% Function implementation
function [D, Q] = im2front_base(I, start_pts, method, ...
    alpha, rho, sigma, der, int, samp, eign )

% we DO NOT accept a variable number of inputs (dealt with in IM2FRONT)
% if nargin<=9,     eign = 'l1';
%     if nargin<=8,     samp = 1;
%         if nargin<=7,     int = 'fast';
%             if nargin<=6,     der = 'fast';
%                 if nargin<=5,     sigma = 1;
%                     if nargin<=4,     rho = 3;
%                         if nargin<=3,     alpha = [1 2];   end
%                     end
%                 end
%             end
%         end
%     end
% end

C = size(I,3);

%%
% build the potential/cost function
if C==1 && any(strcmpi(method,{'pix','pixinv','iso'}))
    P = im2potential(I, method, alpha);
    
elseif C==1 && any(strcmpi(method(1:3),{'grd','iso'}))
    P = im2potential(I, method, alpha, sigma);

else % all other cases
    P = im2potential(I, method, alpha, rho, sigma, der, int, samp, eign);
end

%%
% propagate the front through this potential
[D, Q] = potential2front(P, start_pts);

end % end of im2front_base
