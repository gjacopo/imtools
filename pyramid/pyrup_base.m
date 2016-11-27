%% PYRUP_BASE - Base function for PYRUP.
%
%% Syntax
%    Ip = PYRUP_BASE(I, nlevels, stack, filter, order);
%
%% See also
% Related:
% <PYRUP.html |PYRUP|>, 
% <PYRDOWN_BASE.html |PYRDOWN_BASE|>.
% Called: 
% <EXPAND2D_BASE.html |EXPAND2D_BASE|>.

%% Function implementation
function Ip = pyrup_base(I, nlevels, stack, filter, order)

%%
% setting outputs

if stack,   Ip = cell(nlevels);  end;

%%
% compute the hierarchical representation

Io = I;

for i=1:nlevels
    Io = expand2d_base( Io, filter, order ); 
    
    if stack % store all the levels of the pyramid
        Ip{i} = Io;
    end
end

if ~stack % store only the ultimate levels of the pyramid
    Ip = Io;
end

end % end of pyrup_base
