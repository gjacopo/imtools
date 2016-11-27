%% ENHANCELINE - NOT FINISHED IMPLEMENTING.
%
%% Description
% From [GDS10]:
% "The extraction of potential (dark) network pixels is achieved by initially
% suppressing these networks with a parametric - or rank-min - closing by a
% 2D shaped structuring element (SE), followed by a top-hat by the difference 
% between the closing and the input image, - or top-hat by closing"
% From [SG07]:
% "In mathematical morphology terms, the extraction of dark line pixels can 
% be achieved by initially suppressing these lines with a closing by a disk
% shaped structuring element whose diameter slightly exceeds the width of
% the targeted lines.
% A less restrictive closing consists in considering the intersection (i.e.,
% point-wise minimum) of the closing by all subsets of the disk containing
% a fixed number of pixels, see Fig. 1b. This type of closing is known as a
% parametric closing or rank-min closing owing to its representation in terms
% of a rank filter followed by a point-wise minimum with the input image. 
% This is highlighted by computing the difference between the closing and 
% the input image called top-hat by closing."
%
%% References
% [SG07]  P. Soille and J. Grazzini: "Extraction of river networks from
%      satellite images by combining morphology and hydrology", Proc. CAIP,
%      LNCS, vol. 4673, pp. 636-644, 2007.
%      <http://www.springerlink.com/content/7323nx6774021077/>
%
% [GDS10]  J. Grazzini, S. Dillard and P. Soille: "A new generic method for
%      the semi-automatic extraction of river and road networks in low and
%      mid-resolution satellite images", Proc. of SPIE - Image and Signal
%      Processing for Remote Sensing XVI, vol. 7830, pp. 7830071-10, 2010.
%      <http://spiedigitallibrary.org/proceedings/resource/2/psisdg/7830/1/783007_1>

%% Function implementation
function Z = enhanceline(I, N, k)

% (i) rank-min closing: compute the pointwise minimum of the MM closing by 
% all subsets of a disk (NxN) containing a fixed number (K) of pixels

% create the disk-shaped SE of size N
se = strel('disk',N);
[offsets, heights] = getneighbors(se);
% retrieve the corresponding neighbourhood
s = getnhood(se);
x = size(s,1); 
c = (x+1)/2; offsets = offsets+c;
n =  length(heights); % sum(s(:));
% retrieve the positions of the non null (false) elements in the flat
% disk-shaped SE: the combination of k will be selected among those only
pos = sub2ind([x,x], offsets(:,1), offsets(:,2));

k = min(k,n);
% number of operations: nchoosek(length(pos),k)
v = nchoosek(pos,k);

% compute the closing with the first SE
se = false(size(s));
nhood = false(size(s)); nhood(v(1,:)) = true;
Z = imclose(I,nhood); Z = Z(:);
% compute the minimum over all possible SE of size N with k elements
for i=2:size(v,1)
    nhood = se; nhood(v(i,:)) = true;
    z = imclose(I,nhood);
    Z = min([Z,z(:)], [], 2);
end    
% reshape
Z = reshape(Z,size(I));

% (ii) top-hat by closing: compute the difference between the MM closing and
% the input image
Z = Z - I;

% (iii) threshold top-hat: set to 0 all non-zero responses of the top-hat by 
% closing; reset all other pixels to their value in the original image
Z = I.*(Z==0);


end % end of enhanceline



