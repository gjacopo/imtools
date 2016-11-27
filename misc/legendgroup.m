%% LEGENDGROUP - Provide legends for groups of plot handles. 
%
%% Syntax
%   LEGENDGROUP(handle_plot, legend_plot);
%   [l_h, o_h, p_h, t_s] = LEGENDGROUP(handle_p, legend_p, location, interpreter);
%
%% Inputs
% *|handle_p|* : cell (or array) of graphics objects' handles.
%
% *|legend_p|* : cell of strings; each string in |legend_p{i}| is associated
%     with the group of handles in |handle_p{i}|.
%
% *|location, interpreter|* : see |LEGEND| options; default: |location='Best'|
%     and |interpreter='latex'|.
% 
%% See also  
% Related: 
% <matlab:webpub(whichpath('FIGURE')) |FIGURE|>,
% <matlab:webpub(whichpath('LEGEND')) |LEGEND|>.
% Called: 
% <matlab:webpub(whichpath('LEGEND')) |LEGEND|>,
% <matlab:webpub(whichpath('HGGROUP')) |HGGROUP|>,
% <matlab:webpub(whichpath('SET')) |SET|>,
% <matlab:webpub(whichpath('GET')) |GET|>.

%% Function implementation
function [l_h, o_h, p_h, t_s] = ...
    legendgroup(handle_plot, legend_plot, location, interpreter)

if nargin<4,  interpreter = 'tex';
    if nargin<3,  location = 'Best';  end
end

if ~iscell(handle_plot)
    % typically, in the case we have single handle for each group
    handle_plot = num2cell(handle_plot(:),2);
end

if ~iscell(legend_plot)
    if ~ischar(legend_plot)
        error('legendgroup:inputerror', ...
            'second input must be a string or a cell of strings');
    end
    legend_plot = repmat({legend_plot},numel(handle_plot));
    
elseif ~isequal(numel(handle_plot),numel(legend_plot))
    % one item for each group of objects
    error('legendgroup:inputerror', ...
        'input cells must contain same number of elements');
end

ng = numel(handle_plot);

handleGroup = cell(1,ng);
for i=1:ng
    handleGroup{i} = hggroup;
    set(handle_plot{i},'Parent',handleGroup{i})
    % include these hggroup in the legend
    set(get(get(handleGroup{i},'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
end

handleGroup = cell2mat(handleGroup);
[l_h, o_h, p_h, t_s] = legend(handleGroup, legend_plot, ...
    'Location', location, 'Interpreter', interpreter );

end % end of legendgroup