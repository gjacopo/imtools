%% CANNYEDGEFUSE_BASE - Base function for CANNYEDGEFUSE. 
%
%% Syntax
%     edge = CANNYEDGEFUSE_BASE(I, sig2, sig1, serad, der, hyst);
%
%% See also
% Related: 
% <CANNYEDGEFUSE.html |CANNYEDGEFUSE|>.
% Called: 
% <CANNYEDGE_BASE.html |CANNYEDGE_BASE|>.

%% Function implementation
function edgemap = cannyedgefuse_base(I, sig2, sig1, serad, der, hyst, reduce)

%% 
% setting variables

C = size(I,3);

% ensure to combine two distinct scales
if sig1==sig2,   sig2 = sig1 + 0.5;  end;

% in the case 'emax', the final max of the edge maps of the different
% channels is taken on the output of the combination of the two scales, not
% when computing the edge maps.
if ischar(reduce) && strcmpi(reduce,'eor')
    reducefirst = false;
else
    reducefirst = reduce;
end
 
%%
% compute edges at finest scale: sig1
edgemap = cannyedge_base(I, sig1, der, hyst, reducefirst);

%%
% compute edges at coarsest scale: sig2
edgemap2 = cannyedge_base(I, sig2, der, hyst, reducefirst);

%%
% dilate the edges detected at the coarsest scale
se = strel('disk',serad);  
edgemap2 = imdilate(edgemap2,se);

%%
% compute the final edge map as the set of pixels detected as edges in both
% the finest edge map and the dilated coarsest edge map
edgemap = edgemap & edgemap2;

%%
% special output in multispectral case
if C>1 && ischar(reduce) && strcmpi(reduce,'eor')
    edgemap = max(edgemap, [], 3);
end

end % end of cannyedgefuse_base
