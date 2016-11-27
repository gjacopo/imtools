%% FUZZBOUND_BASE - Base function for FUZZBOUND.
% 
%% Syntax
%     fmap = FUZZBOUND_BASE(fmemb, thmax, thcon, thent);
% 
%% See also
% Related:
% <FUZZBOUND.html |FUZZBOUND|>.

%% Function implementation
function fmap = fuzzbound_base(fmemb, thmax, thcon, thent)

%%
% setting internal variables

[X,Y,nC] = size(fmemb);
ind = 0:X*Y-1;

%%
% main calculation

% the process defining fuzzy boundaries can be done via a slicing
mfm = min(fmemb(:)); % should be 0

fmap = zeros(X,Y,3);

%%
% * an individual location x is selected as belonging to a fuzzy boundary 
% if the value of its maximum fuzzy membership value p_max is less than a
% threshold value r. 

f = permute(reshape(fmemb,[X*Y nC]), [2 1]);
[~,i] = max(f,[],1);

%%
% get the maximum membership value
pmax = f(i+nC*ind);

%%
% store the max value in fmap, 1 component
fmap(:,:,1) = reshape(pmax, [X Y]);
if ~isempty(thmax) && thmax<1 
    fmap(:,:,1) = fmap(:,:,1) <= thmax;
end

%%
% * the confusion index for defining fuzzy boundaries, involves
% two fuzzy membership values for each location : the confusion index is
% evaluated by 1.0 minus the difference between the fuzzy membership values
% of location x belonging to the first most likely and the second most
% likely classes. a threshold T is applied so that location x defines a
% fuzzy boundary if the confusion index is greater than a pre-defined 
% threshold v.

f(i+nC*ind) = mfm-eps; % set to the min
pmax2 = max(f,[],1);

fmap(:,:,2) = 1 - reshape(pmax - pmax2, [X Y]);
if ~isempty(thcon) && thcon<1 
    fmap(:,:,2) = fmap(:,:,2)  <= thcon;
end

%%
% * the measure of entropy, for defining fuzzy boundaries, which uses the 
% complete fuzzy membership values for each location. Large values
% indicate low accuracy in classification, while small values indicate high
% accuracy in classification. Then, it is logical to assert that boundaries 
% usually occur where locations have high degrees of fuzziness, that is,
% big values of entropy.

fmap(:,:,3) = - sum( fmemb .* log (fmemb) / log(2), 3);
if ~isempty(thcon) && thcon>0 
    fmap(:,:,3) = fmap(:,:,3)  > thent;
end

end % end of fuzzbound_base
