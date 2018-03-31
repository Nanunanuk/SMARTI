function A=angles2vector(theta,alpha)
% Rotation matritces create unit vectors from angle of incidence (AOI) & Azimut
% angle

ele = 90 - theta;          % Transform AOI to elevation angle
Rz = [cosd(-alpha) sind(-alpha) 0; -sind(-alpha) cosd(-alpha) 0;0 0 1]; % Rotation about the z-axis
Ry = [cosd(-ele) 0 sind(-ele);0 1 0;-sind(-ele) 0 cos(-ele)]; % Rotation about the y-axis
x_rel = Rz*Ry*[1 0 0]'; % Rotate unit vector [1 0 0] first about the y then about the z axis

A = -x_rel'; % Change direction of the vector