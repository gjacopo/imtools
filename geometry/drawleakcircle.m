%% DRAWLEAKCIRCLE - Create 'leaking' circles.
% 
%% Description
% Create a discrete disk (ie. with points on the regular lattice) with an
% approximated number of holes of approximated size. 
%
%% Syntax
%     DRAWLEAKCIRCLE(R, method, gap, ngap);
%     [coord, A] = DRAWLEAKCIRCLE(R, method, gap, ngap, M, disp);
%
%% Inputs
% *|R|* : radius of the circle.
% 
% *|method|* : string defining the method used for creating holes in the set
%     of lattice points sampled on the disk; it is either:
% 
% * |'grid'| for regular holes distribution,
% * |'rand'| for (pseudo) random holes distribution.
% 
% *|gap|* : size of the holes (in number of pixels/points).
% 
% *|ngap|* : desired number of holes.
% 
% *|M|* : (optional) half-size of the output image A displaying the circle;
%     |M| should be |>(R+1)|; default: |M=R+1|.
% 
% *|disp|* : (optional) flag for drawing (discrete) lines; default: |disp=false|.
%
%% Outputs
% *|coord|* : vector |(N,2)| of the (X,Y)-coordinates of the |N| sampled points.
% 
% *|A|* : output logical image of size |(2*M+1,2*M+1)| with the 'leaking' 
%     circle displayed (set to |true| for points whose coordinates are given 
%     by |coord|).
%
%% See also 
% Related: 
% <DRAWCIRCLE.html |DRAWCIRCLE|>, 
% <DRAWGESTALT.html |DRAWGESTALT|>.
% Called:
% <DRAWLEAKCIRCLE_BASE.html |DRAWLEAKCIRCLE_BASE|>.

%% Function implementation
function varargout = drawleakcircle(R, method, gap, ngap, M, disp)

%%
% setting/checking parameters

error(nargchk(4, 6, nargin, 'struct'));
error(nargoutchk(0, 2, nargout, 'struct'));

% set default
if nargin<6
    if nargout==0,  disp = true;
    else            disp = false;     end 
    if nargin<5 || isempty(M)
        M = R+1;
    elseif M<=R
        warning('drawleakcircle:inputwarning', ...
            'too small output image - enlarged domain');
        M = R+1;
    end
    
end

%% 
% main computation

[coord, A] = drawleakcircle_base(R, method, gap, ngap, M);

%% 
% display

if disp
    figure, plot(coord(:,1), coord(:,2), 'k.');
    center = [M M];
    axis([center(1)-1.1*R center(1)+1.1*R center(2)-1.1*R center(2)+1.1*R])
    axis equal
    % figure, imagesc(A), colormap gray, axis image off
end

if nargout>=1,  varargout{1} = coord;
    if nargout==2,  varargout{2} = A;  end
end

end % end of drawleakcircle
