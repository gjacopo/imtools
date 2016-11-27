%% DRAWGESTALT - Draw simple figures displaying 'Gestalt-like' contours and/or shapes.
% 
%% Description
% Create and/or draw simple binary figures of 'Gestalt-like' contours and/or
% shapes. See [PS06,Pras07].
% 
%% Syntax
%     DRAWGESTALT();
%     [contours, shapes] = DRAWGESTALT();
%     [contours, shapes] = DRAWGESTALT(X, Y);
%     [contours, shapes, h] = DRAWGESTALT(X, Y, disp);
%
%% Inputs
% *|X, Y|* : dimensions of the output images.
% 
% *|disp|* : (optional) flag for displaying the gestalt figures; default:
%     |disp=false|.
%
%% Outputs
% *|contours|* : structure whose fields are matrices representing typical
%     Gestalt-like contours with gaps (holes); those fields are any among:
% 
% * |I|, for a junction of type 1 (endlink1: an interrupt edge connected
%          to a straight contour),
% * |L|, for a junction of type L,
% * |X|, for a multiple junction,
% * |E1|, for incomplete contour,
% * |E2|, for composed incomplete contour,
% * |Y|, for intersection of type Y,
% * |T|, for intersection of type T.
% 
% *|shapes|* : ibid, structure representing typical (opened) shapes, its fields
%     can be any among:
% 
% * |N1|, for single-sided neck (wedge),
% * |N2|, for double-sided neck,
% * |F|, for a tapering (flaring) shape,
% * |Y|, for intersection of type Y,
% * |T|, for intersection of type T,
% * |L|, for intersection of type L.
% 
% *|h|* : when |disp=true|, returns the figure id.
%
%% References
% [PS06]  L. Prasad and A. Skourikhine: "Vectorized image segmentation
%      via trixel agglomeration", Pattern Recognition, 39:501-514, 2006.
%      <http://www.sciencedirect.com/science/article/pii/S0031320305003857>;
%      see also <http://www.springerlink.com/content/yvkn22dxbx6evxnf/>
%
% [Pras07]  L. Prasad: "Rectification of the chordal axis transform skeleton
%      and criteria for shape decomposition", Image and Vision Computing,
%      25:1557-1571, 2007.
%      <http://www.sciencedirect.com/science/article/pii/S026288560600309X>;
%      see also <http://www.springerlink.com/content/hjf8wjvn822cljgt/>
%
%% See also 
% Related: 
% <DRAWCIRCLE.html |DRAWCIRCLE|>.
% <DRAWLEAKCIRCLE.html |DRAWLEAKCIRCLE|>.
% Called:
% <DRAWGESTALT_BASE.html |DRAWGESTALT_BASE|>.

%% Function implementation
function varargout = drawgestalt(X, Y, disp)

%%
% checking/setting parameters

error(nargchk(0, 3, nargin, 'struct'));
error(nargoutchk(0, 3, nargout, 'struct'));

%X=150; Y=200;
if nargin<3,  disp = false;
    if nargin<2,  Y = [];
        if nargin<1,  
            X =150;
            if nargout == 0,  disp = true;  end
        end
    end
end

if nargin==1 && islogical(X)
    disp = X;
    X = 150;
end

gap = round(X/4);

if nargout==0 && ~disp
    error('drawgestalt:inputwarning', 'stupid you... see help and try again');
end

if isempty(Y),  Y = X;  end % default square size

if disp,  
    if ~isempty(ver('images'))
        se = strel('square',3);
        fdisp = @(x) imdilate(x,se);
    else
        fdisp = @(x) x;
    end
    h = figure;  hold on;
else
    if nargout==3,
        warning('drawgestalt:inputwarning', ...
            'no figure displayed - empty ''h'' argument returned');
        h = [];
    end
end

%%
% main calculation

[C, S] = drawgestalt_base(X, Y, gap);

%%
% display

if disp,
    % junction type 1: endlink1
    subplot(4,4,1), imagesc(~fdisp(C.I)), axis image off;
    title('junction type 1 (endlink1)');
    % intersection type L
    subplot(4,4,2), imagesc(~fdisp(C.L)), axis image off;
    title('L junction');
    % intersection type Y
    subplot(4,4,3), imagesc(~fdisp(C.Y)), axis image off;
    title('Y junction');
    % intersection type T
    subplot(4,4,4), imagesc(~fdisp(C.T)), axis image off;
    title('Y junction');
    % multiple junction type 1
    subplot(4,4,5), imagesc(~fdisp(C.X)), axis image off;
    title('multiple junction');
    % incomplete contour
    subplot(4,4,6), imagesc(~fdisp(C.E1)), axis image off;
    title('incomplete');
    % composed incomplete contours
    subplot(4,4,7), imagesc(~fdisp(C.E2)), axis image off;
    title('composed incomplete');
    % single-sided neck
    N = S.N1; N(10:X-10,10) = true;  N(10:round(X/2)+10,Y-10) = true;
    subplot(4,4,9), imagesc(~fdisp(N)), axis image off;
    title('single-sided neck');
    % double-sided neck
    N2 = S.N2;  N2(10:X-10,10) = true;  N2(10:X-10,Y-10) = true;
    subplot(4,4,10), imagesc(~fdisp(N2)), axis image off;
    title('double-sided neck');
    % tapering shape
    F = S.F;  F(10:X-10,10) = true;  F(round(X/2)-10:round(X/2)+10,Y-10) = true;
    subplot(4,4,11), imagesc(~fdisp(F)), axis image off;
    title('tapering shape');
    % intersection type Y
    sY=S.Y;  sY(10:round(X/2)-15,10) = true;  sY(10:round(X/2)-15,Y-10) = true;
    sY(X-10,round(Y/2)-20:round(Y/2)+20) = true;
    subplot(4,4,12), imagesc(~fdisp(sY)), axis image off;
    title('Y-shape');
    % intersection type T
    T = S.T;  T(20:round(X/2),10) = true;  T(20:round(X/2),Y-10) = true;
    T(X-10,round(Y/2)-round(gap/2):round(Y/2)+round(gap/2)) = true;    
    subplot(4,4,13), imagesc(~fdisp(T)), axis image off;
    title('T-shape');
    % intersection type L
    L = S.L;  L(20,20:20+gap) = true;  L(X-10-gap:X-10,round(2*Y/3)) = true;    
    subplot(4,4,14), imagesc(~fdisp(L)),axis image off;
    title('L-shape');
    % main figure's title
    colormap gray, suptitle('''Gestalt-like'' contours and shapes');
    hold off;
end

%%
% final outputs
if nargout>=1,  varargout{1} = C;
    if nargout>=2,  varargout{2} = S;  
        if nargout==3,  varargout{3} = h;  end
    end
end

end % end of drawgestalt