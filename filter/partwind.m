% Creates scaled window of set size.
% function to calculate window centered at the zero segment of outwind size
function windd = partwind(wfunc, insidewind, outwind)
%   function to calculate window centered at the zero segment of outwind size
%   WFUNC - windowing function from the list:
%   'barthannwin' 'bartlett' 'blackman' 'blackmanharris' 'bohmanwin'
%   'gausswin' 'hann' 'nuttallwin' 'parzenwin' 'triang' 'hamming' 'flattopwin'
%   OUTWIND  - length of the OUTPUT total window segment (padded with zero
%   if OUTWIND> INSIDEWIND or truncated INSIDEWIND if OUTWIND < INSIDEWIND
%   INSIDEWIND - length of the windowing function centered in OUTWIND
%
%EXAMPLE:
%       
%       windd1 = partwind('flattopwin', 20, 60);
%       windd2 = partwind('flattopwin', 60);
%       windd3 = partwind('flattopwin', 90, 60);
%       figure; plot(windd1, 'r.-'); hold on; plot(windd2, 'b+-');
%       plot(windd3, 'ko-'); grid; axis tight
%       legend('insidewind<outwind', 'insidewind=outwind', 'insidewind>outwind')

% Other m-files required: none
% Subfunctions: none
% MAT-files required: 
% window.m from Signal Processing Toolboxes

%____________________________________________
% 	Sergei Koptenko, Resonant Medical Inc., 
%           Montreal, Qc., Canada
%	sergei.koptenko@resonantmedical.com 
%   Website: http://www.resonantmedical.com
%____________Apr/22/2005_____________________

if (nargin<3), outwind = insidewind; end;

wnd1 = floor(insidewind/2);       % Get half window size
wnd2 = floor(outwind/2);       % Get half window size
insidewind = 2* wnd1+1; %Force INSIDE window size to be ODD
outwind = 2* wnd2+1; %Force OUTPUT window size to be ODD
eval(['windi =window(@' wfunc ', insidewind);']); % windowing function

%___NEED TO SCALE THIS WINDOW___
if wnd1 < wnd2, % inside window is SMALLER that output
   %disp('wnd1 < wnd2')
    windd = zeros(outwind,1);
    windd(wnd2-wnd1+1:wnd2+wnd1+1) = windi;   
elseif wnd1 > wnd2, % inside window is BIGGER that output
    %disp('wnd1 > wnd2')
    windd= windi(wnd1-wnd2+1:wnd1+wnd2+1);
else     %disp('wnd1 = wnd2')
    windd= windi; % dimensions are equal
end
