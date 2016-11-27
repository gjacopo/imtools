%% TENSCALEDIRFILT_BASE - Base function for TENSCALEDIRFILT.
%
%% Syntax
%     F = TENSCALEDIRFILT_BASE(I, T, S, sig1, sig2, nsig, nthe, nani, aecc, ani);
% 
%% See also
% Related:
% <TENSCALEDIRFILT.html |TENSCALEDIRFILT|>,
% <CONVOLUTION_BASE.html |CONVOLUTION_BASE|>,
% <TENSANIFILT_BASE.html |TENSANIFILT_BASE|>,
% <ADAPTIVEFILT_BASE.html |ADAPTIVEFILT_BASE|>,
% <GEODESICFILT_BASE.html |GEODESICFILT_BASE|>,
% <MDLFILT_BASE.html |MDLFILT_BASE|>.
% Called:
% <ADAPTIVEFILT_BASE.html |ADAPTIVEFILT_BASE|>,
% <../../kernel/html/DIRGAUSSKERNEL.html |DIRGAUSSKERNEL|>.

%% Function implementation
function F = ...
    tenscaledirfilt_base(I, T, S, sig1, sig2, nsig, nthe, nani, aecc, ani)  %#ok

[x y] = size(I(:,:,1));                                         
v = sqrt(x*x + y*y);
m = max(I(:)) - min(I(:));

%%
% possibly build a tensor field
if size(T,4)==2
   [sig1,sig2,~,~,theta] = gstdecomp(T);
    % theta = mod( atan2(e1(:,:,2), e1(:,:,1)), pi );
    
elseif size(T,4)==1
    e1 = perform_vf_normalization(T);
    %  e2 = e1(:,:,2:-1:1); e2(:,:,1) = -e2(:,:,1);
    sig1 = sig1 * ones(size(T,1));
    sig2 = sig2 * ones(size(T,1));
    theta = mod( atan2(e1(:,:,2), e1(:,:,1)), pi );
end

sig1 = v * sig1 / m;
sig2 = v * sig2 / m;

%%
% compute the anisotropy of the tensor field
if ani
    A = (sig1 - sig2) ./ (sig1 + sig2);
    % 
else
    A = sig1 ./ sig2;
end

%%
% compute the sigma product
if isempty(S)
    S = sig1 .* sig2; % square of the sigma deviations
elseif isscalar(S)
    S = S * ones(x,y);
end

if std(A(:))<1e-5,    nani = 1;  end
if std(S(:))<1e-5,    nsig = 1;  end

%% 
% setting the local filters

%%
% compute the ranges of sigma for which different filters need to be built
s = sort(S(:)); a = sort(A(:));
k = round(0.05*length(s));
s_list = linspace(s(k), s(end-k), nsig);
a_list = linspace(a(k), a(end-k), nani);
% quantize the orientation
t_list = linspace(0,pi,nthe+1); t_list(end) = [];

% compute the parameters of the filters

%%
% generation of arrays for 3D function and interpolation
[the_list,sig1_list,ani_list] = ndgrid(t_list,s_list,a_list);
the_list = the_list(:);
ani_list = ani_list(:);

aecc = sqrt(ani_list);
% aecc = (aecc + ani_list(:)) / aecc;

%%
% compute the deviation along short and long axis

% deviation along the short axis
sig2_list = sqrt(sig1_list(:)) ./ aecc;
% deviation along the long axis
sig1_list = sqrt(sig1_list(:)) .* aecc;
m = length(the_list);

%%
% quantize
T = repmat( reshape(the_list,[1 1 m]), [x y 1] );
S1 = repmat( reshape(sig1_list,[1 1 m]), [x y 1] );
S2 = repmat( reshape(sig2_list,[1 1 m]), [x y 1] );
t = repmat( theta, [1 1 m] );
s1 = repmat( sig1, [1 1 m] );
s2 = repmat( sig2, [1 1 m] );

E = abs(T-t) + abs(S1-s1) + abs(S2-s2);
[~,A] = min(E,[],3);

%% 
% performing the adaptive filtering

H =  dirgausskernel(sig1_list, sig2_list, the_list, min(x,y), 'geu');

F = adaptivefilt_base(I, H, A);

end % end of tenscaledirfilt_base

%% Subfunction

%%
% |PERFORM_V_NORMALIZATION| - Renormalize a vector field.
%--------------------------------------------------------------------------
function v = perform_vf_normalization(v)
a = nb_dims(v);
d = sqrt( sum(v.^2,a) );
d(d<eps) = 1;
v = v .* repmat( 1./d, [ones(a-1,1)' size(v,a)] );
end % end of perform_vf_normalization
