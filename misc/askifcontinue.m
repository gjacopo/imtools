%% ASKIFCONTINUE - Utilitiy function to pause a program execution.
%
%% Description
% Pause a program and ask if pursueing it or not.
%
%% Syntax
%         [yes, varargout] = ASKIFCONTINUE(question);
%
%% See also  
% Related: 
% <matlab:webpub(whichpath('PAUSE')) |PAUSE|>.

%% Function implementation
function [yes, varargout] = askifcontinue(question)

if ~exist('question','var') || isempty(question) 
    question = 'continue? y/n [y]: ';
end
reply = input(question, 's');

if isempty(reply) || any(strcmpi(reply,{'y','Y'}))
    yes = true;
elseif any(strcmpi(reply,{'n','N'}))   
    yes = false;
end
for i=1:nargout, varargout{i}=[];  end

end % end of askifcontinue
