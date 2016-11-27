%% FINDIMPORTANTEXTREMA1D_BASE - Find important extrema of 1D signal.
%
%% Description
% Find strict extrema in (1D) time series and  compute their importance
% using the approach of [FG07].
%
%% Syntax
%      [Extrema, Strict, Left, Right] = findimportantextrema1D(X);
%      [Extrema, Strict, Left, Right] = findimportantextrema1D(X, dist);
%
%% Reference
% [FG07]  E. Fink and H.S.Gandhi: "Important extrema of time series", Proc.
%      IEEE International Conference on Systems, Man, and Cybernetics, pp. 
%      366-372, 2007.
%
%% See also
% Related:
% <FINDZEROEXTREMA1D.html |FINDZEROEXTREMA1D|>.
% Called:
% <FINDZEROEXTREMA1D_BASE.html |FINDZEROEXTREMA1D_BASE|>.

%% Function implementation
function [Extrema, Strict, Left, Right] = findimportantextrema1D(X, varargin)

%%
% parsing parameters
error(nargchk(1, 15, nargin, 'struct'));
error(nargoutchk(0, 3, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(X) || nb_dims(X)~=1
    error('findimportantextrema1D:inputerror','a 1D signal is required in input'); 
end

% optional parameters
p = createParser('FINDIMPORTANTEXTREMA1D');   % create an instance of the inputParser class.
p.addOptional('dist', 'abs', @(x) ischar(x) && any(strcmpi(x,{'abs','nabs'})));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%% 
% main computation

[Extrema, Strict, Left, Right] = findimportantextrema1D_base(X, p.dist);

%%
% display

if p.disp
    figure, hold on;
    plot(X,'-kx');  plot(Extrema(:,1), Extrema(:,2), 'ro');
    text(Extrema(:,1), Extrema(:,2), num2str(Strict), 'VerticalAlignment','top','Color','r');
    title('extrema''s strict importance');  hold off;
end

end % end of findimportantextrema1D