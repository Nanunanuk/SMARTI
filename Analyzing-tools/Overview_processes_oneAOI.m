clear all
clc

[FileName,PathName] = ...
 uigetfile(pwd, 'Datensatz auswählen');

cd(PathName);
load(FileName);

A_str='';
for layer=1:size(layer_mat,2)
    A_str=[A_str 'A_' num2str(layer) '(:,2), '];
end
figure

eval(['area(A_1(:,1),[' A_str 'A_AR(:,2), A_bb(:,2)+A_fi(:,2), R_esc(:,2), T_tot(:,2)]);'])
axis([lambda_min lambda_max 0 1])

ylabel('Share of individual pocesses')
xlabel('wavelength \lambda')

legend(strsplit([A_str 'Absorption in AR-coating, ',...
    'Absorption in grid, ',...
    'Reflection, ',...
    'Transmission'],', '),...
    'Location','South')


