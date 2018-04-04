%|************************************************************************
%|SMARTI lambertsch
%|by Nils Reiners
%|This source code is free: you can redistribute it and/or modify it under the terms of the GNU General Public License (GPL) as published by the Free Software Foundation., version 3 of the License. Any redistribution of		 |% %|modified code must also be free under the GNU General Public License.																																							 |%
%|This source code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for |% %|more details. You should have received a copy of the GNU General Public License along with this source code. If not, see <http://www.gnu.org/licenses/>.																		 |%
%|************************************************************************


function V_rand=lambertsch(V)

%% Creating a ray with random direction around the normal vector of the area that was hit
% Input: Row or column vector which represents the normal vector of a surface
% Approach: 1) Creating a random vector around [0,0,1]
%           2) Rotating it by the same sequenz as the vector N need to be
%              rotated to become [0,0,1] but in the reverse direction
% Information on "rotation matrix" can be found on wikipedia

if size(V,1)==3 % Changing to row vector if is column
    V=V';
end

N=V/norm(V);

rand_nr1=randi(1000000)/1000000;
rand_nr2=randi(1000000)/1000000;

%% We create a random vector around the elementary vector [0,0,1]
%  by first rotating it randomly arround the x-axis and the around the z-axis

Rx_rand=[1,0,0
    ;0,cos(-pi/2*rand_nr1),-sin(-pi/2*rand_nr1)
    ;0,sin(-pi/2*rand_nr1),cos(-pi/2*rand_nr1)];

Rz_rand=[cos(2*pi*rand_nr2),-sin(2*pi*rand_nr2),0
    ;sin(2*pi*rand_nr2),cos(2*pi*rand_nr2),0
    ;0,0,1]; % Rotation around the z-axis

%% Reverse rotation as the vector need to be rotated to becom [0,0,1]

if N(2)>=0
    Rz_back=[cos(2*pi-(2*pi+atan(N(1)/N(2)))),-sin(2*pi-(2*pi+atan(N(1)/N(2)))),0
        ;sin(2*pi-(2*pi+atan(N(1)/N(2)))),cos(2*pi-(2*pi+atan(N(1)/N(2)))),0
        ;0,0,1]; % Rotation back around the z-axis
else
    Rz_back=[cos(2*pi-(pi+atan(N(1)/N(2)))),-sin(2*pi-(pi+atan(N(1)/N(2)))),0
        ;sin(2*pi-(pi+atan(N(1)/N(2)))),cos(2*pi-(pi+atan(N(1)/N(2)))),0
        ;0,0,1]; % Rotation back around the z-axis
end

Rx_back=[1,0,0
    ;0,cos(2*pi-(acos(N(3)))),-sin(2*pi-(acos(N(3))))
    ;0,sin(2*pi-(acos(N(3)))),cos(2*pi-(acos(N(3))))]; % Rotation back around the x-axis


if isequal([0,0,1],N)
    V_rand=Rz_rand*Rx_rand*[0;0;1];
else
    V_rand=Rz_back*Rx_back*Rz_rand*Rx_rand*[0;0;1];
end



