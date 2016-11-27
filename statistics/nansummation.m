%% NANSUMMATION - Overload NANSUM.
%
%% Description
% Generic function to overload |NANSUM| in the case the statistics toolbox
% is not available.
%
%% Syntax
%     S = NANSUMMATION(x,dim);
%
%% See also
% Related:
% <NANAVERAGE.html |NANAVERAGE|>.
% Called:
% <matlab:webpub(whichpath('NANSUM')) |NANSUM|>,
% <matlab:webpub(whichpath('SUM')) |SUM|>.
 
%% Function implementation
function S = nansummation(x,dim)

if isempty(ver('stats'))
    x(isnan(x)) = 0;
    S = sum(x,dim);
else
    S = nansum(x,dim);
end
end % end of nansummation