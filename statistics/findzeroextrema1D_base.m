%% FINDZEROEXTREMA1D_BASE - Base function for FINDZEROEXTREMA1D. 
%
%% Syntax
%   [Max, Min, Zero] = FINDZEROEXTREMA1D_BASE(X, x);
%   [Max, Min, Zero] = FINDZEROEXTREMA1D_BASE(X, x, delta);
%   [Max, Min, Zero] = FINDZEROEXTREMA1D_BASE(X, x, met_max, met_min, met_zeros);
%   [Max, Min, Zero] = FINDZEROEXTREMA1D_BASE(X, x, ...
%              met_max, met_min, met_zeros, delta);
%   [Max, Min] = FINDZEROEXTREMA1D_BASE(X, x, met_max, met_min, [], delta);
%   Max = FINDZEROEXTREMA1D_BASE(X, x, met_max, [], [], delta);
%   Min = FINDZEROEXTREMA1D_BASE(X, x, [], met_min, [], delta;
%   Zero = FINDZEROEXTREMA1D_BASE(X, x, [], [], met_zeros, delta);
%
%% Remarks
% * A feature won't be computed only if its method is passed as an empty string:
% eg., no |Max| is computed (and output) when |met_max=[]|; however, a default
% feature is always computed when its method is not passed (e.g. all features
% are output when calling |FINDZEROEXTREMA1D_BASE(X, x, delta))|. 
%
% * |delta| is a |(1,n)| vector, with |n<=3|, of the form |[delta mindelta valabs]|
% (see help |FINDZEROEXTREMA1D| for explanation).
%
% * When no extrema and/or zeros are found, this function returns a row vector
% of null indices (|[0 0]| for |Min| and |Max|, |[0]| for |Zero|).
% This is implemented this way to enable further used of this function besides
% its call in the function |FINDZEROEXTREMA1D| (which, on the contrary, returns
% an empty vector in such case).
%
%% See also
% Related:
% <FINDZEROEXTREMA1D.html |FINDZEROEXTREMA1D|>.
% Called:
% <matlab:webpub(whichpath('PEAKFINDER')) |PEAKFINDER|>,
% <matlab:webpub(whichpath('PEAKDET')) |PEAKDET|>,
% <matlab:webpub(whichpath('FINDEXTREMA')) |FINDEXTREMA|>,
% <matlab:webpub(whichpath('EXTREMA')) |EXTREMA|>,
% <matlab:webpub(whichpath('CROSSING')) |CROSSING|>,
% <matlab:webpub(whichpath('FIND')) |FIND|>.

%% Function implementation
function varargout = ...
    findzeroextrema1D_base(X, x, met_max, met_min, met_zeros, delta )

%%
% check/set internal parameters

if nb_dims(X)>1
    error('findzeroextrema1D_base:inputerror', ...
        '1D signal required in input');
elseif length(X)<3
    for i=1:nargout,  varargout{i} = [0 0];  end
    return
end

% special case
if nargin==3 && isnumeric(met_max)
    delta = met_max;
    met_max = 'peakfinder';

elseif nargin<6, 
    delta = (max(X)-min(X))/4; % see PEAKFINDER
end

% default values for default approach
if nargin<5,  met_zeros = 'crossing'; % by default, we compute it
    if nargin<4,  met_min = 'peakfinder';
        if nargin<3,  met_max = 'peakfinder';
            if nargin<2,  x = []; end
        end
    end
end

if isempty(x),  x = (1:length(X))';  end

if size(X,2)>1,  X = transpose(X);  end
if size(x,2)>1,  x = transpose(x);  end

interpol = false;

% default values when empty options are passed
if length(delta)<3,  valabs = Inf;
    if length(delta)<2,  mindelta = 0;
        if isempty(delta),  delta = 0;  end
    else
        mindelta = delta(2);
    end
else
    valabs = delta(3);
    mindelta = delta(2);
end
delta = delta(1);

%%
% main computation starts here

narg =0;

if any(strcmpi('morph',{met_max,met_min,met_zeros}))
    d = [X(2:end-1) X(3:end) X(1:end-2)];
    I = find(max(abs(diff([d d(:,1)],1,2)),[],2)<delta);
end

if ~isempty(met_max) 
    % no max computation takes place only when met_max is passed as an empty 
    % variable
    narg = narg+1;
    switch met_max
        
         case 'ecke' 
            %% 
            % |ECKE| - Iterative filtering.
            S = X;
            for j = 1:2
                A = diff([X;X(1)]);
                B = flipud(diff(flipud([X(end);X])));
                % min = find(~(((A<0)&(B<=0))|((A<=0)&(B<0))));
                S(~(((A<0)&(B<=0)) | ((A<=0)&(B<0)))) = 0;
                % figure, plot(X, 'x-')
            end
            Max = find(S);
        
        case {'local3','morph'}
            %%
            % |'local3'| - Fast naive method for finding extrema in local
            % (3,1) neighbourhoods.
            % since local maximum and minimum points of a signal have zero
            % derivative, their locations can be estimated from the zero-
            % crossings of |diff(X)|, provided the signal is sampled with
            % sufficiently fine resolution;
            
            % for a coarsely sampled signal, a better estimate is
            A = diff(X(1:end-1)); A(abs(A)<eps) = 0;
            B = diff(X(2:end));  B(abs(B)<eps) = 0;
            % Max = find(sign(A)-sign(B)>0 & abs(A)>delta & abs(B)>delta) + 1;
            Max = find(sign(A)-sign(B)>0 & ...
                min(abs([A B]),[],2)>=delta) + 1;
            
        case 'local5'
            %%
            % |LOCAL5| - Method for finding extrema in local (5x1) neighbourhoods.
            % typically, each signal value is compared to its 4 neighbours 
            % (2 left, 2 right) to check if it is a local extremum.
            % examples of accepted configurations ('o' represents the tested
            % value, the '_' represents a locally constant signal, while '/'
            % represents a locally increasing signal and '\' represents a
            % locally decreasing signal:
            %
            %      o       /o\   _/o\      o__    o_     o
            %    _/ \_    /   \      \    /      /  \   / \_
            %                            /      /      /
            % examples of rejected configurations
            %
            %     _o_       _     /
            %    /   \    o/    o/   /\o/\  /\o       o      o_/
            %           _/     /               \/    / \/   /
            %                 /                     /      /
            A = diff(X(2:end-2));  A(abs(A)<eps) = 0;
            B = diff(X(3:end-1));  B(abs(B)<eps) = 0;
            C = diff(X(1:end-3));  C(abs(C)<eps) = 0;
            D = diff(X(4:end));  D(abs(D)<eps) = 0;
            Max = find(sign(A)-sign(B)>0 & sign(C)-sign(D)>0 & ...
                min(abs([A B]),[],2)>=delta) + 2;
           
           
        case 'peakfinder'
            %%
            % |PEAKFINDER| - Noise tolerant fast peak finding algorithm.
            % |PEAKFINDER| quickly finds local peaks or valleys (local extrema)
            % in a noisy vector using a user-defined magnitude threshold to
            % determine if each peak is significantly larger (or smaller) than
            % the data around it. The problem with the strictly derivative
            % based peak finding algorithms is that if the signal is noisy
            % many spurious peaks are found. However, more complex methods
            % often take much longer for large data sets, require a large
            % amount of user interaction, and still give highly variable
            % results. This function attempts to use the alternating nature
            % of the derivatives along with the user defined threshold to
            % identify local maxima or minima in a vector quickly and robustly.

            [Max, vMax] = peakfinder(X, delta, 1);
            % [peakLoc, peakMag] = peakfinder(x0,thresh,extrema) returns the
            % indices of the local maxima as well as the magnitudes of those
            % maxima.
            %   INPUTS:
            %       x0 - A real vector from the maxima will be found (required)
            %       thresh - The amount above surrounding data for a peak to
            %           be identified (default = (max(x0)-min(x0))/4). Larger
            %           values mean the algorithm is more selective in finding
            %           peaks.
            %       extrema - 1 if maxima are desired, -1 if minima are desired
            %           (default = maxima, 1)
            %   OUTPUTS:
            %       peakLoc - The indicies of the identified peaks in x0
            %       peakMag - The magnitude of the identified peaks
            
            
        case 'peakdet'
            %%
            % |PEAKDET| - Detect peaks in a vector.
            % Using the well-known zero-derivate method. Due to the noise,
            % which is always there in real-life signals, accidental zero-
            % crossings of the first derivate occur, yielding false detections.
            % The typical solution is to smooth the curve with some low-pass
            % filter, usually killing the original signal at the same time.
            % The result is usually that the algorithm goes horribly wrong
            % where it's so obvious to the eye.
            % The trick here is to realize, that a peak is the highest point
            % between "valleys". What makes a peak is the fact that there are
            % lower points around it. This strategy is adopted by |PEAKDET|:
            % look for the highest point, around which there are points lower
            % by X on both sides.
            
            [Max, Min] = peakdet(X, delta); % Max and Min are already complete
           % [maxtab, mintab] = peakdet(v, delta, x) finds the local maxima
            % and minima ("peaks").
            %   INPUTS:
            %       V - input vector.
            %       delta - a point is considered a maximum peak if it has
            %           the maximal value, and was preceded (to the left) by
            %           a value lower by delta: we require a difference of 
            %           at least delta between a peak and its surrounding in
            %           order to declare it as a peak.
            %       x - the indices in MAXTAB and MINTAB (see below) will be
            %           replaced with the corresponding x-values.
            %   OUTPUTS:
            %       MAXTAB, MINTAB - consist of two columns; column 1 contains
            %           indices in V, and column 2 the found values.
            
            
        case 'fpeak'
            %%
            % |FPEAK| - Find peak value of data.
            % |FPEAK| reports the minimum data points and requires more inputs 
            % than the |PEAKFINDER| method for instance.

            varargout{1} = fpeak(x, X, delta);
            % peak = fpeak(x, y, s, Range)
            %   INPUTS:
            %       x, y - input signal coordinates and values
            %       s - sensitivity of the function.
            %       Range - peak value's range
            varargout{2} = [];
            
        case 'findextrema'
            %%
            % |FINDEXTREMA| - find indices of local extrema and zero-crossings.
            [Max, Min, Zero] = findextrema(X);
             % [IMAX,IMIN,ICRS] = FINDEXTREMA(X) returns the indices of local
            % maxima in IMAX, minima in IMIN and zero-crossing in ICRS for 
            % input vector X            
           
        case 'extrema'
            %%
            % |EXTREMA| - Gets the global extrema points from a time series.
            % |EXTREMA| reports many maxima peaks
            [vMax,Max,vMin,Min] = extrema(X);
            % [XMAX,IMAX,XMIN,IMIN] = EXTREMA(X) returns the global minima 
            % and maxima points.
            %   INPUTS:
            %       X - input vector, possibly containing NaN values (ignored)
            %   OUTPUTS:
            %       XMAX - maxima points in descending order
            %       IMAX - indexes of the XMAX
            %       XMIN - minima points in descending order
            %       IMIN - indexes of the XMIN
                       
    end
    
    if strcmpi(met_max,'morph')
        S = X;
        ero = min(d,[],2); ero(ero<0) = 0;
        S(I+1) = ero(I);
        A = diff(S(1:end-1)); B = diff(S(2:end));
        M = find(sign(A)-sign(B)>0) + 1;
        Max = Max(ismember(Max,M));
    end
    
    if ~isnan(mindelta) && mindelta>0
        if ~strcmpi(met_max,'local3')
            A = diff(X(1:end-1)); B = diff(X(2:end));
        end
        md = find(max(abs([A B]),[],2)>mindelta)+1;
        Max(~ismember(Max,md)) = [];
        if exist('vMax','var') && ~isempty(vMax),  
            vMax(ismember(Max,md)) = [];
        end
    end
   
    if ~isnan(valabs) && valabs<Inf && valabs>0,  
        mv = find(abs(X)<=valabs); 
        Max(ismember(Max,mv)) = [];
        if exist('vMax','var') && ~isempty(vMax),  
            vMax(ismember(Max,mv)) = [];
        end
    end
 
    if ~strcmpi(met_max,'fpeak')
        % retrieve the coordinates of the extrema in the original system and
        % determine the extrema values, if not done already
        if size(Max,2)<2,
            if ~exist('vMax','var') || isempty(vMax),  vMax = X(Max);  end
            varargout{narg} = [x(Max) vMax];
        else
            varargout{narg} = Max;
            varargout{narg}(:,1) = x(Max(:,1));
        end
    end
    
    if isempty(varargout{narg}), varargout{narg} = [0 0];  end
end

if ~isempty(met_min)
    % no min computation takes place only when met_min is passed as an empty
    % variable
    narg = narg+1;
    switch met_min
        
        case 'ecke' % nothing output
            S = X;
            for j = 1:2
                A = diff([X;X(1)]);
                B = flipud(diff(flipud([X(end);X])));
                S(~(((A>0)&(B>=0)) | ((A>=0)&(B>0)))) = 0;
            end
            Min = find(S);
            
        case {'local3','morph'}
            if isempty(met_max) || ~strcmpi(met_max,'local3')
                A = diff(X(1:end-1));  A(abs(A)<eps) = 0;
                B = diff(X(2:end));  B(abs(B)<eps) = 0;
                % else: A and B have been already previously computed
            end
            % Min = find(sign(A)-sign(B)<0 & abs(A)>delta & abs(B)>delta) + 1;
            Min = find(sign(A)-sign(B)<0 & ...
                min(abs([A B]),[],2)>=delta & ...
                max(abs([A B]),[],2)>mindelta) + 1;
            
        case 'local5'
            if isempty(met_max) || ~strcmpi(met_max,'local5')
                A = diff(X(2:end-2));  A(abs(A)<eps) = 0;
                B = diff(X(3:end-1));  B(abs(B)<eps) = 0;
                C = diff(X(1:end-3));  C(abs(C)<eps) = 0;
                D = diff(X(4:end));  D(abs(D)<eps) = 0;
            end
            Min = find(sign(A)-sign(B)<0 & sign(C)-sign(D)<0 & ...
                min(abs([A B]),[],2)>=delta) + 2;
            
        case 'peakfinder'
            [Min, vMin] = peakfinder(X,delta,-1);
            
        case 'peakdet'
            if isempty(met_max) || ~strcmpi(met_max,'peakdet')
                [~, Min] = peakdet(X, delta);
                % else: we already computed a Min
            end
            
        case 'fpeak'
            if isempty(met_max) || ~strcmpi(met_max,'fpeak')
                varargout{narg} = fpeak(x, X, delta);
            end
            
        case 'findextrema'
            if isempty(met_max) || ~strcmpi(met_max,'findextrema')
                [~, Min, Zero] = findextrema(X);
                % else: we already computed a Min
            end
            
        case 'extrema'
            if isempty(met_max) || ~strcmpi(met_max,'extrema')
                [~, ~, vMin, Min] = extrema(X);
                % else: we already computed Min, vMin
            end
            
    end
    
    if strcmpi(met_max,'morph')
        S = -X; d = -d;
        ero = min(d,[],2); ero(ero<0) = 0;
        S(I+1) = ero(I);                
        A = diff(S(1:end-1)); B = diff(S(2:end));
        m = find(sign(A)-sign(B)>0) + 1;
        Min = Min(ismember(Min,m));
    end
    
    if ~isnan(mindelta) && mindelta>0
        if ~strcmpi(met_min,'local3')
            A = diff(X(1:end-1)); B = diff(X(2:end));
        end
        md = find(max(abs([A B]),[],2)>mindelta)+1;
        Min(~ismember(Max,md)) = [];
        if exist('vMin','var') && ~isempty(vMin),  
            vMin(ismember(Min,md)) = [];
        end
    end
    
    if ~isnan(valabs) && valabs<Inf && valabs>0,  
        if ~exist('mv','var'),  mv = find(abs(X)<=valabs);  end
        Min(ismember(Min,mv)) = [];
        if exist('vMin','var') && ~isempty(vMin),  
            vMin(ismember(Min,mv)) = []; 
        end
    end
    
    if ~strcmpi(met_min,'fpeak')
        if ~isempty(Min)
            if size(Min,2)<2,
                if ~exist('vMin','var') || isempty(vMin),  vMin = X(Min);  end
                varargout{narg} = [x(Min) vMin];
            else
                varargout{narg} = Min;
                varargout{narg}(:,1) = x(Min(:,1));
            end
        else
            varargout{narg} = [];
        end
    end
    if isempty(varargout{narg}), varargout{narg} = [0 0];  end
end

if strcmpi(met_zeros,'morph')
    if any(strcmpi('morph',{met_max,met_min}))
        X = S;
    else
        ero = min(d,[],2); ero(ero<0) = 0;
        X(I+1) = ero(I);
    end
end

if ~isempty(met_zeros)
    narg = narg+1;
    switch met_zeros
        
        case 'naive'
            %%
            % Naive method
            % to obtain  the indices where signal x crosses zero
            i0 =  find(diff(sign(X)));
            % the kth zero-crossing lies between x(i(k))  and x(i(k)+1)
            % linear interpolation can be used for subsample estimates of
            % zero-crossings locations
            if interpol     % linear interpolation
                varargout{narg} = i0 - X(i0)./(X(i0+1) - X(i0));       %#ok
            else
                varargout{narg} = i0;
            end
            
        case {'local3','morph'}
            %%
            % |'local3'| - Extraction of zero crossing by comparing a value
            % with its direct neighbours.
            A = diff(X(1:end-1));  B = diff(X(2:end));
            % same thing, but then B=-B:
            %  A = diff(X);  A = A(1:end-1);
            %  B = flipud(diff(flipud(X)));  B = B(2:end);
            % thus, we then have to check A.*B<=0 in the following
            varargout{narg} = ...
                find(A.*B>=0 & X(1:end-2).*X(3:end)<=0 & ... % zero-crossing
                min(abs([A B]),[],2)>delta) ...              % non constant
                + 1;
            % & abs(X(2:end-1))<= min([abs(X(1:end-2)) abs(X(3:end))],[],2)
            
        case 'crossing'
            %%
            % |CROSSING| - Find the crossings of a given level of a signal
            if interpol
                [ind,varargout{narg}] = crossing(X);                   %#ok
            else
                varargout{narg} = crossing(X);
            end
            % [ind,t0] = CROSSING(S,t,level,par)
            %   INPUTS:
            %       S - input signal.
            %       t - interpolating time (possibly [] for no interpolation).
            %       level - crossings at this level will be returned instead
            %           of the zero crossings.
            %       par - optional string {'none'|'linear'}: with interpolation
            %	        turned off (par = 'none') this function always returns
            %	        the value left of the zero (the data point that is
            %           nearest to the zero AND smaller than the zero
            %           crossing).
            %   OUTPUTS:
            %       ind - index vector so that the signal S crosses zero at
            %           ind or between ind and ind+1
            %       t0 - time vector t0 of the zero crossings of the signal
            %           S; the crossing times are linearly interpolated
            %           between the given times t
            
        case 'findextrema'
            %%
            % |FINDEXTREMA| - find indices of local extrema and zero-crossings.
            if any(strcmpi({met_max; met_min},{'findextrema';'findextrema'}))
                % Zero has been computed already
                varargout{narg} = Zero;
            else
                [~, ~, varargout{narg}] = findextrema(X); % see above
            end
            
    end
    
    if ~isnan(valabs) && valabs~=Inf,
        if ~exist('mv','var'),   mv = find(abs(X)<=valabs);  end
        varargout{narg}(~ismember(varargout{narg},mv)) = [];
    end
    
    if ~isnan(mindelta) && mindelta>0
        if ~strcmpi(met_zeros,{'local3','morph'})
            A = diff(X(1:end-1)); B = diff(X(2:end));
        end
        md = find(max(abs([A B]),[],2)>=mindelta)+1;
        varargout{narg}(~ismember(varargout{narg},md)) = [];
    end
   
    % retrieve the coordinates in the original system
    varargout{narg} = x(varargout{narg});
    if isempty(varargout{narg}), varargout{narg} = 0;  end
end

if narg<nargout  % too many variables have been required in output
    varargout(narg+1:nargout) = cell(nargout-narg,1);
    if narg==0
        error('findzeroextrema1D_base:errorinput', ...
            ['at least one non empty string must be entered as a method - ' ...
            'see variables met_min, met_max and met_zeros']);
    end
end

end % end of findzeroextrema1D_base