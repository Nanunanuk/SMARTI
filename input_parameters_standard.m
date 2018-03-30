%|************************************************************************
%|SMARTI SETTING INPUT-PARAMETERS
%|by Nils Reiners
%|This source code is free: you can redistribute it and/or modify it under the terms of the GNU General Public License (GPL) as published by the Free Software Foundation., version 3 of the License. Any redistribution of		 |% %|modified code must also be free under the GNU General Public License.																																							 |%
%|This source code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for |% %|more details. You should have received a copy of the GNU General Public License along with this source code. If not, see <http://www.gnu.org/licenses/>.																		 |%
%|************************************************************************


%% The following 6 parameters can be used for batch simulation [Value1 Value2 ...]
batch_1=[0];
batch_2=[0];
batch_3=[75];
batch_4=[0];
batch_5=[54.74];    
batch_6=[100]; 

%% DIRECTION OF INCIDENT LIGHT
alpha = batch_1;
theta = batch_2;

rdm_ray_dir = 0; % Diffuse sky (1 to activate)

%% Setting spectral resolution
lambda_min=300;
lambda_max=1200;
lambda_step=30;

%% Setting number of rays (first range < 960nm; second for > 960nm)
nr_of_rays_1=1000;
nr_of_rays_2=1000;

%% Loading material parameter. Use names from excel sheet in Data-Folder => {'GeneralName', 'SpecificName', [Source]}
% More materials can be loaded than the actually used ones
Mat(1,:)={'Air' 'Air' '[Pal85a]'};
Mat(2,:)={'Si' 'Crystalline 300 K' '[Gre08]'};
Mat(3,:)={'AlSi' 'Eutectic' '[Vog15]'};
Mat(4,:)={'Ag','Pure','[Pal85a]'};
Mat(5,:)={'SiNx','PECVD','[Bak11]'};
Mat(6,:)={'Glass','Low-Fe Starphire','[McI09b]'};
Mat(7,:)={'EVA','UV transmissive','[Vog16]'};

load_data(Mat)
Reflector=0.5; % Reflector can also be chosen as material below
%% Determination of the materials (Use the 'General' discription of the material in the Mat file)
layer_mat={Air Glass EVA Si AlSi};

%% Setting material for coating on top of corresponding layer (if no coating is applied write 'noAR')
AR_mat={'noAR' 'noAR' SiNx 'noAR' 'noAR'};
d_AR={0 0 batch_3 batch_3 0};                    %nm

%% Choose the geometry of each layer (Use the name of the corresponding geometry function starting with 'geometry')
S_geom={...
    'geometry_cube_pyramide'...
    'geometry_cube_pyramide'...
    'geometry_cube_pyramide'...
    'geometry_cube_pyramide'...
    'geometry_cube_pyramide'}
w=10; %[µm] Width of the unit cell
p_bot={w/2*tand(batch_4) w/2*tand(batch_5) 0 0 0}; % Depth of texture at bottom of layer

%% SUBSTRATE POSITION
Substrate_pos=4;                    % Substrate position is needed for the generation profile

d_layer={10 2000 30 200 10};       %µm

% Lambertian scattering on bottom of the layer => Insert value between 0 and 1
lambert={0 0 0.4 0.4 0};

%% GRID PARAMETERS
grid.mat={Ag Ag};

pf=2150; %Pitch between fingers [µm]
wf=180; %Width of the fingers [µm]

pbb=46; % [mm]
wbb=2.48;% [mm] The busbar share can be set to zero by inserting a low value such as 0.0001

grid.pos={0 0 1 0 0}; %On top of the layer. No grid if all are 0.

%% DISPLAY 3D-PLOT OF UNIT CELL (WITH OR WITHOUT RAYS)
plotting_geom=0;
plotting_rays=0;
xth_ray=20;

%% Detection of the transmission angles after the refraction at the adjacent medium
detect_bot=[0 0 1 0 0]; % Example: If you want to detect the angles the rays are entering the active lyer (e.g. silicon) then set the 'detect_bot' value in the layer tha is on top of it to 1
detect_top=[0 0 0 0 0]; 

%% WRITING GENERATION PROFILES (as .gen-Files for PC1D)
gen_create=0;