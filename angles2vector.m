%|************************************************************************
%|SMARTI angles2vector
%|by Nils Reiners
%|This source code is free: you can redistribute it and/or modify it under the terms of the GNU General Public License (GPL) as published by the Free Software Foundation., version 3 of the License. Any redistribution of		 |% %|modified code must also be free under the GNU General Public License.																																							 |%
%|This source code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for |% %|more details. You should have received a copy of the GNU General Public License along with this source code. If not, see <http://www.gnu.org/licenses/>.																		 |%
%|************************************************************************

function A=angles2vector(theta,alpha)
% Rotation matritces create unit vectors from angle of incidence (AOI) & Azimut
% angle

ele = 90 - theta;          % Transform AOI to elevation angle
Rz = [cosd(-alpha) sind(-alpha) 0; -sind(-alpha) cosd(-alpha) 0;0 0 1]; % Rotation about the z-axis
Ry = [cosd(-ele) 0 sind(-ele);0 1 0;-sind(-ele) 0 cos(-ele)]; % Rotation about the y-axis
x_rel = Rz*Ry*[1 0 0]'; % Rotate unit vector [1 0 0] first about the y then about the z axis

A = -x_rel'; % Change direction of the vector