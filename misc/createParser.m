%% CREATEPARSER - Create an instance of the inputParser class.
%
%% Syntax
%       p = CREATEPARSER(fname);
%
%% Input 
%   fname : string containing the name of the function where the parser is
%     created.
%
%% Output
%   p : an instance of the inputParser class.
%
%% See also  
% Related: 
% <GETVARPARSER.html |GETVARPARSER|>,
% <matlab:webpub(whichpath('INPUTPARSER')) |INPUTPARSER|>.

%% Function implementation
function p = createParser(fname)

if ~ischar(fname)
    error('createParser:inputargument',...
        'input argument must a string containing the function name');
end

p = inputParser;   % create an instance of the inputParser class.
% single structure 'parse' passing all the (optional) parameters: always
% present in the parser
p.addParamValue('parse', [], @isstruct);
% additional default options
p.addParamValue('verb', false, @islogical);
p.addParamValue('debug',false, @islogical);
p.addParamValue('disp', false, @islogical);
% name of the function
p.FunctionName = fname;
p.CaseSensitive = true;

end % end of createParser