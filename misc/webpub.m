%% WEBPUB - Creation of web formatted functions' documentation.
%
%% Syntax
%     WEBPUB(name, update, outputDir, showCode, evalCode);
%
%% Inputs
% See |PUBLISH| options.
%
% *|name|* : a string or a cell of strings specifying either:
%   
% * a directory name whose all contained files will be published, 
% * the name of a file to publish,
% * the names of several files to publish.
%
% *|update|* : an (optional) logical flag forcing, when |true| the update of
%     the published document (html) when it already exists (and have been
%     found in the user path); default: |update=false|.
%
% *|outputDir|* : (optional) string specifying the folder where to publish 
%     the output documentation and its associated image files:
%
% * |''| (default): output is placed in the html subfolder of the function 
%          folder,
% * full path: output is placed in the specified folder.
%
% *|showCode|* : (optional) logical value that specifies whether the code is
%     to be included in the published documentation; default: |showCode=true|.
%
% *|evalCode|* : (optional) logical value specifying whether to run the code
%     that is published; default: |evalCode=false|.
%
%% Output
% Create an html document in the subdirectory _html_ of the directory where
% the file(s) to be published was(were) found (_html_ is created if it does
% not already exist) and display the created document.
%
%% See also links
% Related:
% <matlab:web(whichpath('HELP')) |HELP|>,
% <matlab:web(whichpath('DOC')) |DOC|>.
% Called:
% <matlab:web(whichpath('WEB')) |WEB|>,
% <matlab:web(whichpath('PUBLISH')) |PUBLISH|>.

%% Function implementation
function webpub(name, update, outputDir, showCode, evalCode, verbose)

%%
% check/set variables

if nargin<6,  verbose = false;
    if nargin<5,  evalCode = false;
        if nargin<4,  showCode = true;
            if nargin<3,  outputDir = '';  
                if nargin<2,  update = false; end
            end
        end
    end
end

%%
% define the list of names of the file to publish
if ischar(name)
    if isdir(name) % retrieve all the M-files in the directory
        if verbose,  disp('M-files in the directory to be published');  end
        listing = dir([name '/*.m']);
        if isempty(listing)
            error('webdoc:errorinput', ...
                '!!!directory not found - set full (or relative) path!!!');
        end
        filename = {listing(:).name};
    else % transform the string in a cell of string
        if verbose,  disp('single file to be published');  end
        filename = {name};
    end    
elseif iscell(name)
    if verbose,  disp('list of files to be published');  end
    filename = name;
    
else
    error('webpub:errorinput', ...
        'Argument to webpub must be a string or a cell of strings');
end
%%
% * create the publishing options

opt.outputDir = outputDir;
opt.showCode = showCode;
opt.evalCode = evalCode;
opt.format = 'html';

%%
% loop over the file names for publishing all of them

nfile = length(filename);
if verbose && nfile>1
    disp([num2str(nfile) ' files found to be published']);  
elseif nfile==0,  
        error('webdoc:errorinput', '!!!no file found!!!');
end

for i=1:nfile
    %%
    % * check that the file exists (possible in cases where no directory was entered)
    if ~exist(filename{i},'file')
        error('webdoc:errorinput', ...
            ['!!!file/function ' filename{i} ' not found!!!']);
    end
    
    %%
    % * get rid of the '.m' extension (could use STRFIND instead...)
    if strcmpi(filename{i}(end-1:end),'.m'),
        filename{i} = filename{i}(1:end-2);
    end
    
    %%
    % * edit after possible prior publishing
    if ~update && exist([filename{i} '.html'],'file')
        % the file has been already published and it is in a known path
        tmp = [filename{i} '.html']; % simply edit
    else
        if verbose,  disp(['publishing file ' filename{i} ' ...']);  end
        tmp = publish(filename{i},opt); % publish
    end
end

%%
% web display the last published file in the case a directory was passed as
% argument
if verbose,  disp(['editing file ' tmp]);  end
web(tmp);

end % end of webpub
