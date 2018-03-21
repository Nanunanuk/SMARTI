function A=angles2vector(theta,alpha)
% Drehmatritzen erzeugen Vektor aus Elevation & Azimut

ele = 90 - theta;          % Zenith in Elevation umrechnen
Rz = [cosd(-alpha) sind(-alpha) 0; -sind(-alpha) cosd(-alpha) 0;0 0 1]; % Drehung um z-Achse
Ry = [cosd(-ele) 0 sind(-ele);0 1 0;-sind(-ele) 0 cos(-ele)];
x_rel = Rz*Ry*[1 0 0]';

A = -x_rel';