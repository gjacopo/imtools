%% PYRDOWN_BASE - Base function for PYRDOWN.
%
%% Syntax
%    Ip = PYRDOWN_BASE(I, nlevels, stack, filter, order);
%
%% See also
% Related:
% <PYRDOWN.html |PYRDOWN|>, 
% <PYRUP_BASE.html |PYRUP_BASE|>.
% Called: 
% <REDUCE2D_BASE.html |REDUCE2D_BASE|>.

%% Function implementation
function Ip = pyrdown_base(I, nlevels, stack, filter, order)

%%
% setting outputs

if stack,    Ip = cell(nlevels);  end;

%% 
% compute the hierarchical representation

Io = I;

for i=1:nlevels
    Io = reduce2d_base( Io, filter, order ); 
    
    if stack % store all the levels of the pyramid
        Ip{i} = Io;
    end
end

if ~stack % store only the ultimate levels of the pyramid
    Ip = Io;
end

end % end of pyrdown_base
