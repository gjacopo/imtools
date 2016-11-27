%% NANAVERAGE - Overload NANMEAN.
%
%% Description
% Generic function to overload |NANMEAN| in the case the statistics toolbox
% is not available.
%
%% Syntax
%     S = NANAVERAGE(x,dim);
%
%% See also
% Related:
% <NANSUMMATION.html |NANSUMMATION|>.
% Called:
% <matlab:webpub(whichpath('NANSUM')) |NANMEAN|>,
% <matlab:webpub(whichpath('SUM')) |SUM|>.

%% Function implementation
function M = nanaverage(x,dim)

if isempty(ver('stats'))
nans = isnan(x);
    x(nans) = 0;
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    M = sum(x,dim) ./ n;
else
    M = nanmean(x,dim);
end

end % end of nanaverage