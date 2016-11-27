%% GETVARPARSER - Return an instance  of the inputParser class.
%
%% Syntax
%       param = GETVARPARSER(p);
%
%% Input 
%   p : an instance of the inputParser class.
% 
%% Output
%   param : full list of parameters contained in the input Parser structure.
%
%% See also  
% Related: 
% <CREATEPARSER.html |CREATEPARSER|>,
% <matlab:webpub(whichpath('INPUTPARSER')) |INPUTPARSER|>.

%% Function implementation
function param = getvarParser(p)

% get the variables entered with a 'parse' structure (in p.Results.parse)
param = p.Results.parse;

% fill unexisting fields with other entered values (in p.Results)
fields = fieldnames(p.Results);
for ip=1:numel(fields)
    if ~isfield(param,fields{ip})
        param.(fields{ip}) = p.Results.(fields{ip});
    end
end

% param.FunctionName = p.FunctionName

% old approach
% p = p.Results;
% fields = fieldnames(p)
% for ip=1:numel(fields)
%     eval([genvarname(fields{ip}) '= p.(fields{ip});']);
% end
% % possibly overwrite variables when a 'param' variable has been passed,
% % while unpassed default values are kept
% if ~isempty(param)                                                   
%     fields = fieldnames(param);
%     for ip=1:numel(fields)
%         eval([genvarname(fields{ip}) '= param.(fields{ip});']);
%     end
% end
% p=[]; param=[];                                                      

% special case (predefined verbose and debug fields)
if param.debug,  param.verb = true;  param.disp = true;  end

if param.debug,
    disp(['in ' p.FunctionName ':']);
    disp(param);    
end

end % end of getvarParser
    