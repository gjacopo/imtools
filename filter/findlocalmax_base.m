%% FINDLOCALMAX_BASE - Base function for FINDLOCALMAX.
%
%% Syntax
%     [max_local,row,col] = findlocalmax_base(I, radius, method );  
%
%% See also
% Related:
% <FINDLOCALMAX.html |FINDLOCALMAX|>,
% <FINDLOCALEXTREMA_BASE.html |FINDLOCALEXTREMA_BASE|>.
% Called: 
% <matlab:webpub(whichpath('ORDFILT2')) |ORDFILT2|>,
% <matlab:webpub(whichpath('IMDILATE')) |IMDILATE|>,
% <matlab:webpub(whichpath('SORT')) |SORT|>.

%% Function implementation
function [max_local,row,col] = findlocalmax_base(I, radius, method)

%% 
% dealing with multispectral images

[M,N,C] = size(I);
if C>1
    max_local = zeros(size(I)); 
    row = cell(C); col = cell(C);
    for c=1:C
        [max_local(:,:,c), row{c}, col{c}] = ...
            findlocalmax_base(I(:,:,c), radius, method);
    end
    return;
end

%% 
% find the local maxima

% define the neighbourhood
mask  = fspecial('disk',radius)>0;
    
if any(strcmp(method,{'dil','filt'}))
    switch method
        case 'dil' % find local maxima by dilation (fast) /!\ non unique /!\
            I2 = imdilate(I,mask);
            index = I==I2;
            
        case 'filt'  % find unique local maxima using filtering (fast)
            nb    = sum(mask(:));
            highest          = ordfilt2(I, nb, mask);
            second_highest   = ordfilt2(I, nb-1, mask);
            index            = highest==I & highest~=second_highest;
    end
    
    max_local        = zeros(size(I));
    max_local(index) = I(index);
    [row,col]        = find(index==1);

elseif strcmp(method, 'screen') % find unique local maxima (slow)
    max_local   = zeros(M,N);
    I_enlarge = zeros(M+2*radius,N+2*radius);
    I_mask    = zeros(M+2*radius,N+2*radius);
    I_enlarge( (1:M)+radius , (1:N)+radius ) = I;
    I_mask(    (1:M)+radius , (1:N)+radius ) = 1;
    row = zeros(M*N,1);
    col = zeros(M*N,1);
    index = 0;
    for l = 1:M
        for c = 1:N
            I_ref = I(l,c);
            neigh_I  = I_enlarge(l:l+2*radius,c:c+2*radius);
            neigh_mask = I_mask(   l:l+2*radius,c:c+2*radius).*mask;
            neigh_sort = sort(neigh_I(neigh_mask==1));
            if I_ref==neigh_sort(end) && I_ref>neigh_sort(end-1)
                index          = index+1;
                row(index,1)   = l;
                col(index,1)   = c;
                max_local(l,c) = I_ref;
            end
        end
    end
    row(index+1:end,:) = [];
    col(index+1:end,:) = [];
end

end % end of findlocalmax_base
