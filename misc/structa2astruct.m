%% STRUCTA2ASTRUCT - Convert a structure of arrays to an array of structures.
%
%% Syntax
%          A = STRUCTA2ASTRUCT(S);
%          [A, fields] = STRUCTA2ASTRUCT(S, destroy, fields, ...);
% 
%% Remark
% See discussion at
% http://blogs.mathworks.com/pick/2008/04/22/matlab-basics-array-of-structures-vs-structures-of-arrays/
%
%% See also  
% Related: 
% <matlab:webpub(whichpath('STRUCT')) |STRUCT|>,
% <matlab:webpub(whichpath('STRUCT2CELL')) |STRUCT2CELL|>.

%% Function implementation
function [A,fields] = structa2astruct(S,destroy,varargin)

error(nargoutchk(1, 2, nargout, 'struct'));
if ~isstruct(S)
    error('structa2astruct:inputerror','structure required in input') 
end

if isempty(destroy), destroy=false; end

fields = fieldnames(S);

for ip=1:numel(fields)
    
    if nargin>1 && ~any(strcmpi(fields{ip},varargin))
        continue;
        % else: we agree to fill A with all fields found in S
    end
    
    if iscell(S.(fields{ip})) % already a cell array
        cell_field = S.(fields{ip});
        if ~isempty(cell_field)
            % assign the field across the structure array,
            [A(1:length(cell_field)).(fields{ip})] = cell_field{:};
        end
    
    elseif isstruct(S.(fields{ip})) % a structure: recursive call
        warning('structa2astruct:method','unsolved conversion of struct')
        struct_field = stra2astr(S.(fields{ip}));                      %#ok        
        
    else  % create cells (concatenating the columns)
        cell_field = num2cell(S.(fields{ip}),2);
        if ~isempty(cell_field)
            % assign the field across the structure array,
            [A(1:length(cell_field)).(fields{ip})] = cell_field{:};
        end
    end
        
    % possibly clean up
    if destroy,  S.(fields{ip}) = [];  end
    
end        

end % end of structa2astruct