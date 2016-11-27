%% ASF_BASE - Base function for ASF.
%
%% Syntax
%     F = ASF_BASE(I, n, shape, op, s1, s2);
%     F = ASF_BASE(I, n, shape, op, s1, s2);
%
%% See also
% Related:
% <ASF.html |ASF|>,
% <GRANULOMETRY_BASE.html |GRANULOMETRY_BASE|>,
% <MORPHPROFILE_BASE.html |MORPHPROFILE_BASE|>.
% Called:
% <matlab:web(whichpath('IMOPEN')) |IMOPEN|>,
% <matlab:web(whichpath('IMCLOSE')) |IMCLOSE|>,
% <IMROPEN.html |IMROPEN|>,
% <IMRCLOSE.html |IMRCLOSE|>,
% <FLATSTREL.html |FLATSTREL|>.

%% Function implementation
function F = asf_base(I, n, shape, op, s1, s2)

%%
% preparing internal operators and parameters

if ischar(shape)
    if nargin>=5 && isnumeric(s1)
        if nargin<6 || ~isnumeric(s2),  s2 = 0;  end
    else
        s1 = 0;  s2 = 0;
    end
    
elseif isnumeric(shape)
    tmp = shape;
    if nargin>=5 && ischar(s1),  shape = s1;
        if nargin<6 || ~isnumeric(s2),  s2 = 0;  end
    else
        shape = []; s2 = 0;
    end
    s1 = tmp;
end

%%
% main computation

if strcmp(op(1),'r')
    fopen = @imropen;
    fclose = @imrclose;
    op = op(2:end);
else
    fopen = @imopen;
    fclose = @imclose;
end

switch op
    case 'oc'
        first_op = fclose;
        second_op = fclose;
        third_op = @(x,y)x;
        
    case 'co'
        first_op = fclose;
        second_op = fopen;
        third_op = @(x,y)x;
        
    case 'oco'
        first_op = fopen;
        second_op = fclose;
        third_op = fopen;
        
    case 'coc'
        first_op = fclose;
        second_op = fopen;
        third_op = fclose;
        
end

F = I;

for i=1:n
    se = flatstrel(shape, s1, s2);
    
    F = third_op(second_op(first_op(F,se), se), se);
    
    switch shape
        case {'line','square','rectangle'}
            s1 = s1 + 2;                                               
        case {'diamond','disk','octagon','periodicline'}
            s1 = s1 + 1;                                               
    end
end

end % end of asf_base
