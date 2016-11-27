%% MORPHPROFILE_BASE - Base function for MORPHPROFILE.
%
%% Syntax
%     DMP = MORPHPROFILE_BASE( I, op, se );
%     DMP = MORPHPROFILE_BASE( I, op, se );
%     [DMP, Phi, MP] = MORPHPROFILE_BASE( I, op, se, s1[, s2] );
%
%% See also
% Related:
% <MORPHPROFILE.html |MORPHPROFILE|>,
% <GRANULOMETRY_BASE.html |GRANULOMETRY_BASE|>,
% <ASF_BASE.html |ASF_BASE|>.
% Called: 
% <IMRECONSTRUCTBY.html |IMRECONSTRUCTBY|>,
% <FLATSTREL.html |FLATSTREL|>.
 
%% Function implementation
function [DMP, Phi, varargout] = morphprofile_base( I, op, se, s1, s2 )

[X,Y] = size(I(:,:,1));   

n = min(length(op),3);

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

R0 = I;

if any(strcmp(op,{'roc','oc','ed'}))
    DMP = cell(1,2*N); % note: in principle, MP{0}=I
else
    DMP = cell(1,N);
end

if strcmp(op(1:n),'roc'),
    R00 = I;
    if nargout>2,
        switch op
            case 'roc',      varargout{1} = cell(2*N,1);
            case 'rocmax',   varargout{1} = cell(N,1);
        end
    end
end

%% 
% compute the (derivative of the) morphological profile

for i=1:N
    
    if strcmp(op(1),'r')
        switch op
            case {'roc','rocmax'}
                R = imreconstructby_base(I, 'ro', se{i});  % Eq.(1) of [PB]
                R2 = imreconstructby_base(I, 'rc', se{i}); % Eq.(2)               
            case {'rc','rclose','ro','ropen'}
                R = imreconstructby_base(I, op, se{i});
        end
        
    else     
        switch op(1)
            case {'o','open'}
                R = imopen(I, se{i});
            case {'c','close'}
                R = imclose(I, se{i});
            case {'e','erode'}
                R = imerode(I, se{i});
            case {'d','dilate'}
                R = imdilate(I, se{i});
        end
        if length(op)>=2
            if strcmp(op(1:2),'oc'),  R2 = imclose(I, se{i});
            elseif strcmp(op(1:2),'ed'),  R2 = imdilate(I, se{i});
            end
        end
        
    end
    
    %     if der
    DMP{i} = abs(R - R0); % Eq.(3)
    if nargout>2,  varargout{1}{N+i} = R;  end
    if any(strcmp(op,{'roc','oc','ed'}))
        DMP{2*N-i+1} = abs(R2 - R00); % Eq.(4)
    elseif any(strcmp(op,{'rocmax','ocmax','edmax'}))
        DMP{i} = max(cat(3, DMP{i}, R2 - R00), [], 3);
    end
    if nargout>2 && any(strcmp(op,{'roc','oc','ed','rocmax','ocmax','edmax'})),
        varargout{1}{2*N-i+1} = R2;
    end
    
    %     else
    %         MP{i} = R;
    %         if any(strcmp(op,{'roc','oc','ed'})),
    %             MP{2*N-i+1} = R2;
    %         end
    %     end
    
    R0 = R;
    if any(strcmp(op,{'roc','oc','ed','rocmax','ocmax','edmax'}))
        R00 = R2;
    end

end

%% 
% compute the morphological multi-scale characteristics

Phi = zeros(size(I));

T = cell2mat(DMP);

for c = 1:size(T,3)
    [~, Phi(:,:,c)] = max(reshape(T(:,:,c),[X,Y,length(DMP)]), [], 3);
end

end % end of morphprofile_base
