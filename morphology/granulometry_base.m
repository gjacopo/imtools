%% GRANULOMETRY_BASE - Base function for GRANULOMETRY.
%
%% Syntax
%     G = GRANULOMETRY_BASE(I, op, se, s1, s2);
% 
%% Remark
% This is a non-recursive implementation!
% 
%% See also
% Related:
% <GRANULOMETRY.html |GRANULOMETRY|>,
% <MORPHPROFILE_BASE.html |MORPHPROFILE_BASE|>,
% <ASF_BASE.html |ASF_BASE|>.
% Called: 
% <matlab:web(whichpath('IMOPEN')) |IMOPEN|>,
% <matlab:web(whichpath('IMCLOSE')) |IMCLOSE|>,
% <IMROPEN.html |IMROPEN|>,
% <IMRCLOSE.html |IMRCLOSE|>,
% <FLATSTREL.html |FLATSTREL|>.

%% Function implementation
function G = granulometry_base(I, op, se, s1, s2)

%%
% preparing internal operators and parameters

if nargin>=3 && iscell(se) && strcmp(class(se{1}),'strel')
    N = size(se,1);
    % ignore all other parameters
    
elseif nargin>=3
    if ischar(se)
        shape = se;
        if nargin>=4 && isnumeric(s1)
            N = length(s1);
            if nargin<5 || ~isnumeric(s2),  s2 = cell(N,1);  end
        else
            N = 10;  s1 = cell(N,1);  s2 = cell(N,1);
        end
        
    elseif isnumeric(se)
        N = length(se);
        if nargin>=4 && ischar(s1),  shape = s1;
            if nargin<5 || ~isnumeric(s2),  s2 = cell(N,1);  end
        else
            shape = []; s2 = cell(N,1); 
        end
        s1 = se;
    end   
    
    if ~iscell(s1),  s1 = num2cell(s1);  end
    if ~iscell(s2),  s2 = num2cell(s2);  end
    
    se = cell(N,1);
    for i=1:N
        se{i} = flatstrel(shape, s1{i}, s2{i});
    end
    
end

% if ischar(shape)
%     if nargin>=5 && isnumeric(s1)
%         if nargin<6 || ~isnumeric(s2),  s2 = 0;  end
%     else
%         s1 = 0;  s2 = 0;
%     end
%     
% elseif isnumeric(shape)
%     tmp = shape;
%     if nargin>=5 && ischar(s1),  shape = s1;
%         if nargin<6 || ~isnumeric(s2),  s2 = 0;  end
%     else
%         shape = []; s2 = 0;
%     end
%     s1 = tmp;
% end

G = cell(1,N);

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
    case 'o'
        first_op = fopen;
        second_op = @(x,y)x;
        third_op = second_op;
        
    case 'c'
        first_op = fclose;
        second_op = @(x,y)x;
        third_op = second_op;

    case 'oc'
        first_op = fopen;
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
        
    case 'e'
        first_op = @imerode;
        second_op = @(x,y)x;
        third_op = second_op;
        
    case 'd'
        first_op = @imdilate;
        second_op = @(x,y)x;
        third_op = second_op;
end

% G{1} = I;

for i=1:N
    % se = flatstrel(shape, size, s2);
    
    % se{i}
    G{i} = third_op(second_op(first_op(I,se{i}), se{i}), se{i});
    % G{i} = G{i} - first_op(third_op(second_op(I,se{i}), se{i}), se{i});
    % G{i} = G{i} - fopen(I,se{i});
    %figure, imagesc(rescale(G{i})), colormap gray
    
    % switch shape
    %     case {'line','square','rectangle'}
    %         size = size + 2*s1;                                      %#ok
    %     case {'diamond','disk','octagon','periodicline'}
    %         size = size + s1;                                        %#ok
    % end
    % if i<N,  G{i+1} = G{i};  end;
end

end % end of granulometry_base
