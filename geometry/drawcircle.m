%% DRAWCIRCLE - Create and draw simple circle.
% 
%% Description
% Create and draw simple circle defined by its center, its radius and the
% number of sampled points on it.
%
%% Syntax
%     DRAWCIRCLE(center, radius, NOP);
%     [coord, A] = DRAWCIRCLE(center, radius, NOP, disp, style);
%
%% Inputs
% *|center, radius, NOP|* : resp. center, radius and number of sampled points
%     of the circle to be represented.
% 
% *|disp|* : (optional) flag for drawing (discrete) lines; default: |disp=false|.
% 
% *|style|* : (optional) string defining the style for display (see options for
%     |PLOT|); default: |style = 'b-'|.
% 
%% Outputs
% *|coord|* : |(NOP,2)| matrix storing the (X,Y)-coordinates of the |NOP| grid
%     points used to represent the circle. 
% 
% *|A|* : appproximate representation of the output circle in a logical image.
%
%% See also 
% Related: 
% <DRAWLEAKCIRCLE.html |DRAWLEAKCIRCLE|>, 
% <DRAWGESTALT.html |DRAWGESTALT|>.
% Called:
% <DRAWCIRCLE_BASE.html |DRAWCIRCLE_BASE|>.

%% Function implementation
function varargout = drawcircle(center, radius, NOP, disp, style)

%%
% checking/setting parameters
error(nargchk(2, 4, nargin, 'struct'));
error(nargoutchk(0, 2, nargout, 'struct'));

% set default
if nargin<5
    style = 'b-';
    if nargin<4
        if nargout==0,  disp = true;
        else            disp = false;     end
        if nargin<3,  NOP = 50;  end
    end
end

%% 
% main computation

[coord, A] = drawcircle_base(center, radius, NOP);

if disp,  figure, axis square, plot(coord(:,1), coord(:,2), style);  end

if nargout>=1,  varargout{1} = coord;
    if nargout==2,  varargout{2} = A;  end
end

end % end of drawcircle
