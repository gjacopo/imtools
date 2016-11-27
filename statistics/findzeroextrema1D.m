%% FINDZEROEXTREMA1D - Find extremas and/or zero-crossings of 1D signals.
%
%% Description
% Wrapping function for popular Matlab functions calculating extremas and/or
% zero-crossings of 1D signals. 
%
%% Syntax
%   [Max, Min] = FINDZEROEXTREMA1D(S);
%   [Max, Min, Zero] = FINDZEROEXTREMA1D(S, x, met_extrema, met_zeros, D);
%
%% Inputs
% *|S|* : input 1D signal.
%
% *|x|* : (optional) domain coordinates of the input signal; default: |x=[]|, 
%     ie. the natural coordinates |(1:length(S))| are considered.
%
% *|met_extrema|* : (optional) string defining the function used to compute
%     the extrema of the signal; it is any of the following strings:
%
% * |'fpeak'|, |'peakfinder'|, |'peakdet'|, |'extrema'| or |'findextrema'|
%          when a function with same name needs to be called (see acknowledgment
%          below), 
% * |'local3'| and |'local5'| where a direct methods based on local (3- or 5-)
%          neighbourhoods comparison are implemented, 
% * |'morph'| where a pseudo morphological approach adopted to prior filter
%          the input signal, 
% * |'ecke'| where a two-steps direct approach is used to also prior filter
%          the input (see function |ECKE| in |SHAPEPARTS|); 
%
% default: |met_extrema='local5'|.
%
% *|met_zeros|* : (optional) string defining the function used to compute the 
%     zero-crossings of the signal; it is any among those strings: |'local3'|,
%     |'findextrema'| or |'crossing'|; see below for the corresponding 
%     functions; default: |met_zeros='local3'|. 
%
% *|D|* : (optional) vector storing the variables |[delta mindelta absval]|
%     (in this order) used simultaneously for finding extrema and zero-crossings
%     when (see above) the method |met_extrema| is any among: |'local3'|, 
%     |'morph'|, |'local5'|, |'peakfinder'|, |'peakdet'|, |'fpeak'|,  and/or
%     when the method |met_zeros| is any among |'local3'|, |'morph'|; a point
%     is considered a peak if:
%
% * it is a local extrema, 
% * it is surrounded by neighbour values both low(high)er by at least delta
%          (absolute difference),
% * at least one of those values is low(high)er by mindelta (thus, we should
%          have mindelta>delta), and
% * its absolute value is larger than absval;
%
% a point is considered a zero-crossing if:
%
% * its neighbour values have opposite signs,
% * its absolute value is lower than absval, 
% * its absolute difference with both its neighbour values is larger than
%          delta, and
% * at least one of those differences is larger than mindelta; 
%
% default: |D=[]|, ie. no restriction.
%
%% Outputs
% *|Min, Max|* : matrix storing consisting of two columns; the first column
%     contains the coordinates of the extrema positions, and column 2 contains
%     the found values; in the case |x=[]| or |x=(1:length(X))|, the first
%     column is in fact the indices of the extrema positions in the input 
%     signal; note moreover that in the case |met_extrema='fpeak'|, the minima
%     and maxima are stored together in the matrix |Max| (|Min| is therefore 
%     empty).
% 
% *|Zero|* : coordinates (or indices) of the zero-crossings.
%
%% Acknowledgment
% The external functions called by |FINDZEROEXTREMA1D| are:
%
% * |PEAKFINDER| (noise tolerant fast peak finding algorithm) -
%    from: http://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder
%    see also: http://blogs.mathworks.com/pick/2009/11/13/this-peaks-my-interest/
% * |PEAKDET| (find the local maxima and minima in some noisy signal) -
%    from: http://www.billauer.co.il/peakdet.html
%    see also: http://blogs.mathworks.com/pick/2008/05/09/finding-local-extrema/
% * |FPEAK| -
%    from: http://www.mathworks.com/matlabcentral/fileexchange/4242
% * |FINDEXTREMA| (find indices of local extrema and zero-crossings)
%    from: http://www.mathworks.com/matlabcentral/fileexchange/24306-findextrema
% * |EXTREMA| -
%    from: http://www.mathworks.com/matlabcentral/fileexchange/12275
%    see also: http://blogs.mathworks.com/pick/2008/05/09/finding-local-extrema
% * |CROSSING| (find the crossings at a given level of a signal) -
%    from: http://www.mathworks.com/matlabcentral/fileexchange/2432
%
%% See also
% Related:
% <matlab:webpub(whichpath('PEAKFINDER')) |PEAKFINDER|>,
% <matlab:webpub(whichpath('PEAKDET')) |PEAKDET|>,
% <matlab:webpub(whichpath('FINDEXTREMA')) |FINDEXTREMA|>,
% <matlab:webpub(whichpath('EXTREMA')) |EXTREMA|>,
% <matlab:webpub(whichpath('CROSSING')) |CROSSING|>,
% <matlab:webpub(whichpath('FIND')) |FIND|>.
% Called:
% <FINDZEROEXTREMA1D_BASE.html |FINDZEROEXTREMA1D_BASE|>.

%% Function implementation
function [Max, Min, Zero] = findzeroextrema1D(S, varargin)

%%
% parsing parameters
error(nargchk(1, 15, nargin, 'struct'));
error(nargoutchk(0, 3, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(S) || nb_dims(S)~=1
    error('findzeroextrema1D:inputerror','a 1D signal is required in input'); 
end

% optional parameters
p = createParser('FINDZEROEXTREMA1D');   % create an instance of the inputParser class.
p.addOptional('x', [], @(x) isempty(x) || (isnumeric(x) && nb_dims(x)==1));
p.addOptional('met_extrema','local5', @(x) isempty(x) || (ischar(x)&& ...
    any(strcmpi(x,{'ecke','local3','local5','findextrema','extrema', ...
    'peakfinder','fpeak','peakdet','morph'}))));
p.addOptional('met_zeros','local3', @(x) isempty(x) || (ischar(x)&& ...
    any(strcmpi(x,{'local3','naive','findextrema','crossing','morph'}))));
p.addOptional('delta', [], @(x) isempty(x) || ...
    (length(x)<=3 && isfloat(x) && all(x>=0)));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%%
% setting internal variables

if ~isempty(p.met_extrema) && ~any(strcmpi(p.met_extrema,{'local3','local5','ecke','morph'})) && ...
    ~exist(p.met_extrema,'file')
    error('findzeroextrema1D:unknown', ...
        ['unknown method/function ' met_extrema]);
elseif ~isempty(p.met_zeros) && ~any(strcmpi(p.met_zeros,{'local3','naive','morph'})) && ...
    ~exist(p.met_zeros,'file')
    error('findzeroextrema1D:unknown', ...
        ['unknown method/function ' met_zeros]);
end

if ~isempty(p.x) && ~isequal(size(p.x),size(S))
        error('findzeroextrema1D:inputerror', ...
            'signal and signal''s coordinates must have same dimension'); 
elseif isempty(p.x)
    p.x = (1:length(S));
end

if isempty(p.delta),  p.delta = [0.5 0.1 0.1 0];  end

%% 
% main computation

if nargout<3 && (nargout~=0 || ~p.disp)
    [Max, Min] = findzeroextrema1D_base(S, p.x, ...
        p.met_extrema, p.met_extrema, [], p.delta);
else
    [Max, Min, Zero] = findzeroextrema1D_base(S, p.x, ...
        p.met_extrema, p.met_extrema, p.met_zeros, p.delta);
end

Zero = nonzeros(Zero);

%%
% display

if p.disp
    figure, plot(p.x, S, 'x-'), hold on, grid on;
    plot(Max(:,1), Max(:,2), 'ro');
    if ~isempty(Min),  plot(Min(:,1), Min(:,2), 'gx');  end
    title(['extrema found using ' p.met_extrema]);
    if exist('Zero','var') && ~isempty(Zero)
        figure, plot(p.x, S, 'x-'), hold on, grid on;
        plot(Zero(:,1), S(Zero(:,1)), 'ro');
        title(['zero-crossings found using ' p.met_zeros]);     
    end
end

end % end of findzeroextrema1D