%% COMPRESSRANGE - Compress the entries of a matrix.
% 
%% Description
% Compress (reduce) the range of a matrix of positive (or null) values.
%
%% Syntax
%     I = COMPRESSRANGE(I);
%
%% Input/output
% *|I|* : image with integer values to compress.
%
%% See also
% Called:
% <matlab:webpub(whichpath('CUMSUM')) |CUMSUM|>,
% <matlab:webpub(whichpath('DIFF')) |DIFF|>,
% <matlab:webpub(whichpath('RESHAPE')) |RESHAPE|>.

%% Function implementation
function I = compressrange(I)
 
sz = size(I);
I = I(:);

if ~(isinteger(I) || isequal(I,round(I))) || ~all(I>=0)
    error('compressrange:errorinput', ...
        'only input image with positive integer values are dealt with')
end

values = unique(I);  % sorted list of unique values
if values(1)~=0, values = [0; values];  end

%%
% check if we really have to do something
if length(values)==max(I(:))+1,  return;  end

%%
% this function ensures that the values displayed in the input image are set
% in sequential order 0 values in the original image are left with the label 0
A = cumsum(diff(values)-1);
values = values(2:end);
A = values - A;
% do the relabelling: compress
for l=1:length(values),    I(I==values(l)) = A(l);  end

%%
% reshape
I = reshape(I,sz);

end % end of compressrange