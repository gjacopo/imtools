function [Zones,Omega] = localorientzone(Theta,nz)
% LOCALORIENTZONE - Define the clockwise zone of an image of orientation
%
% Inputs:
%   Theta : image of local orientation in the range [-pi,pi].
%   nz : number of zones used to quantize the orientation.
%
% Outputs:
%   Zones : image providing the representation of the orientation Theta in
%      the range [1,nz].
%   Omega: ratio factor (nz*Theta) / (2*pi).

% We adopt the following index representation for the orientation (see 
% also function LOCAL3X3KERNEL):
% nz=8:
%          ---------------         -------------
%          | NW | N | NE |         | 4 | 3 | 2 |
%          ---------------         -------------
%          | W  |   | E  |    =>   | 5 |   | 1 |
%          ---------------         -------------
%          | SW | S | SE |         | 6 | 7 | 8 |
%          ---------------         -------------
% nz=16:
%          -----------------------         -----------------
%          |  \ NNW | N | NNE /  |         | \ 6 | 5 | 4 / |
%          |   NW   |   |   NE   |         |  7  |   |  3  |
%          | WNW \  |   |  / ENE |         | 8 \ |   | / 2 |
%          -----------------------         -----------------
%          |  W     |   |     E  |    =>   | 9   |   |   1 |
%          -----------------------         -----------------
%          | WSW /  |   |  \ ESE |         | 10/ |   | \16 |
%              SW   |   |   SE   |         |  11 |   | 15  |
%          |  / SSW | S | SSE \  |         | / 12| 13|14 \ |
%          -----------------------         -----------------

% divide the range [-pi,pi] circle
ztheta = linspace(0,2*pi,nz+1);
ztheta(ztheta>pi) = ztheta(ztheta>pi) - 2*pi;
% in the case nz=8, we have:
%      ztheta = [0 0.7854 1.5708 2.3562 3.1416 -2.3562 -1.5708 -0.7854 0]

% define the zones
Zones = ceil(nz*Theta/(2*pi) + nz*(Theta<=0));

%unique(Zones)
% note : theta=0 is sent to zone=nzones
ZTheta = ztheta(Zones+1);

R = nz * abs(ZTheta-Theta) / (2*pi);
Omega = R; % 1 - R; % check
% S = 1 - (1-sqrt(2.)) * R;
end