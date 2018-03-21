%|************************************************************************
%|SMARTI
%|by Nils Reiners
%|This source code is free: you can redistribute it and/or modify it under the terms of the GNU General Public License (GPL) as published by the Free Software Foundation., version 3 of the License. Any redistribution of		 |% %|modified code must also be free under the GNU General Public License.																																							 |%
%|This source code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for |% %|more details. You should have received a copy of the GNU General Public License along with this source code. If not, see <http://www.gnu.org/licenses/>.																		 |%
%|************************************************************************

% Detector für beliebige laser erweitern: testen
% Gen profile off
% R_ext und R_esc trennen
% Test auf beliebige Anzahl von layern erweitern
% Andere Funktionen beschriften


clc
clear all

% Load input parameters
run('input_parameters.m');


% Choose location where the results shell be saved
save_file_loc = uigetdir(pwd,'Select path to store output data');

batch=0;
%% 
for batch_1=batch_1
    for batch_2=batch_2
        for batch_3=batch_3
            for batch_4=batch_4
                for batch_5=batch_5
                    for batch_6=batch_6
                        batch=batch+1;
                        
                        %% Creating folder for results
                        formatOut = 'yyyy_mm_dd_hh_MM';
                        folderName=[datestr(now,formatOut) '_gen_batch-' num2str(batch) '_angle-' num2str(theta)];
                        mkdir(save_file_loc,folderName)
                        Gen_folder=[save_file_loc '\' folderName];
                        
                        %% Creating array of wavelengths
                        lambda=[lambda_min:lambda_step:lambda_max];
                       
                        %% FINGERS AND BUSBARS
                        Mbb=wbb/(wbb+pbb);
                        Mf=pbb*wf/((wbb+pbb)*(wf+pf));
                        Msi=pbb*pf/((wbb+pbb)*(wf+pf));
                        
                        grid.share={Mbb Mf};
                                                                       
                        %% ERROR CHECKING
                        if [p_bot{1}]>[d_layer{1}] || [p_bot{2}]>[d_layer{2}]|| [p_bot{3}]>[d_layer{3}]|| [p_bot{4}]>[d_layer{4}]
                            errordlg('Hight of geom. feature is higher then layer thickness!')                                % Fehlermeldung
                            return
                        end
                        
                        if isequal(layer_mat{Substrate_pos},Si)~=1
                            errordlg('The substrate needs to be silicon!')                                % Fehlermeldung
                            return
                        end
                        %% Creating Cell array including all geometry data
                        layer_data = struct('layer_position',num2cell([1:size(layer_mat,2)]),'AR_mat',AR_mat,'d_AR',d_AR,'layer_mat',layer_mat,'d_layer',d_layer,'p_bot',p_bot,'lambert',lambert);
                        
                        S=0;
                        for S_loop=1:size(layer_mat,2)
                            eval(['S=' S_geom{S_loop} '(' num2str(S_loop) ', S, w, layer_data, plotting_geom);']);
                            S_all{S_loop}=S;
                        end
                                                                   
                        %% STARTING LOOP FOR EACH WAVELENGTH
                        for lambda_batch=[lambda_min:lambda_step:lambda_max] % Loop for every wavelength
                            
                            if lambda_batch<960
                                nr_of_rays=nr_of_rays_1;
                            else
                                nr_of_rays=nr_of_rays_2;
                            end
                            
                            %% Initializing ray matrix
                            I_top=init_rays(rdm_ray_dir,theta,alpha,nr_of_rays,w,lambda_batch,lambda_batch,lambda_step);%lambda_min,lambda_max,lambda_step);
                            I_bot=struct('nr',[],'wavelength',[],'start_loc',[],'dir_vector',[],'ray_col',[]);
                            R_int=struct('nr',[],'wavelength',[],'start_loc',[],'dir_vector',[]);
                            T_tot_struct=struct('nr',[],'wavelength',[],'start_loc',[],'dir_vector',[],'ray_col',[]);
                            error_tot_struct=struct('nr',[],'wavelength',[],'start_loc',[],'dir_vector',[],'ray_col',[],'layer',[],'nr_of_bounces',[]);
                            interrupt_tot_struct=struct('nr',[],'wavelength',[],'start_loc',[],'dir_vector',[],'ray_col',[],'layer',[]);                       
                            A_tot_struct=struct('nr',[],'wavelength',[],'layer_position',[],'abs_loc',[]);                           
                            
                            %% STARTING LOOP FOR EACH RAY
                            for j=1:100; % Loop for every ray
                                
                                
                                %% Layer loop:
                                for layer_pos=1:size(layer_mat,2) % Loop for every layer
                                    
                                    if j>1
                                        if layer_pos~=size(layer_mat,2)
                                            I_bot=eval(['T_top' num2str(layer_pos+1)]);
                                        else
                                            I_bot=struct('nr',[],'wavelength',[],'start_loc',[],'dir_vector',[],'ray_col',[]);
                                        end
                                    end
                                    
                                    [T_top, T_bot, A, error, interrupt]=layer(batch,theta,j,w,layer_data,layer_pos,grid, S_all{layer_pos}, I_top, I_bot, plotting_rays,xth_ray);
                                    I_top=T_bot;
                                    
                                    
                                    if layer_pos==1
                                        % Reflection data
                                        if isempty(R_int(1).nr)==1
                                            R_int=T_top;
                                        elseif isempty(T_top(1).nr)==1
                                            R_int=R_int;
                                        else
                                            R_int=[R_int T_top];
                                        end
                                    else
                                        eval(['T_top' num2str(layer_pos) '=T_top;'])
                                    end
                                    
                                    if layer_pos==size(layer_mat,2)
                                        if isempty(T_tot_struct(1).nr)==1
                                            T_tot_struct=T_bot;
                                        elseif isempty(T_bot(1).nr)==1
                                            T_tot_struct=T_tot_struct;
                                        else
                                            T_tot_struct=[T_tot_struct T_bot];
                                        end
                                        
                                        I_top=struct('nr',[],'wavelength',[],'start_loc',[],'dir_vector',[],'ray_col',[]);
                                    end
                                    
                                    % Absorption data
                                    if isempty(A_tot_struct(1).nr)==1
                                        A_tot_struct=A;
                                    elseif isempty(A(1).nr)==1
                                        A_tot_struct=A_tot_struct;
                                    else
                                        A_tot_struct=[A_tot_struct A];
                                    end
                                    
                                    
                                    %% Inserting errors from layer function in error_tot_struct
                                    if isempty(error(1).nr)~=1&&isempty(error_tot_struct(1).nr)==1
                                        error_tot_struct=error;
                                    elseif isempty(error(1).nr)~=1&&isempty(error_tot_struct(1).nr)~=1
                                        error_tot_struct=[error_tot_struct error];
                                    end
                                    
                                    %% Inserting interrupts from layer function in interrupt_tot_struct
                                    if isempty(interrupt(1).nr)~=1&&isempty(interrupt_tot_struct(1).nr)==1
                                        interrupt_tot_struct=interrupt;
                                    elseif isempty(interrupt(1).nr)~=1&&isempty(interrupt_tot_struct(1).nr)~=1
                                        interrupt_tot_struct=[interrupt_tot_struct interrupt];
                                    end
                                    
                                    %% The detector function is calculating the angles of incidence
                                    if detect_bot(1)==1&&j==1
                                        angle_stats_bot_1=detector(T_bot);
                                        title('Layer 1 in')
                                    end
                                    if detect_top(1)==1&&j==1;
                                        angle_stats_bot_2=detector(T_top);
                                        title('Layer 1 out')
                                    end
                                end

                            end
                            
                            %% ********* PROCESSING RESULTS PER WAVELENGTH ***********
                            row=(lambda_batch-lambda_min)/lambda_step+1;
                            R_esc(row,1:2)=[lambda_batch length([R_int.wavelength])/nr_of_rays];
                            
                            A_fi(row,1:2)=[lambda_batch length([A_tot_struct([A_tot_struct.layer_position]==-1).wavelength])/nr_of_rays];
                            A_bb(row,1:2)=[lambda_batch length([A_tot_struct([A_tot_struct.layer_position]==-2).wavelength])/nr_of_rays];
                            A_AR(row,1:2)=[lambda_batch length([A_tot_struct([A_tot_struct.layer_position]==0).wavelength])/nr_of_rays];
                            
                            for abs_layer=1:size(layer_mat,2)
                                eval(['A_' num2str(abs_layer) '(row,1:2)=[lambda_batch length([A_tot_struct([A_tot_struct.layer_position]==' num2str(abs_layer) ').wavelength])/nr_of_rays];']);
                            end
                                                           
                            T_tot(row,1:2)=[lambda_batch length([T_tot_struct.wavelength])/nr_of_rays];
                            
                            error_tot(row,1:2)=[lambda_batch length([error_tot_struct.wavelength])/nr_of_rays];
                            interrupt_tot(row,1:2)=[lambda_batch length([interrupt_tot_struct.wavelength])/nr_of_rays];
                            
                            %% Creating generation profile
                            start=0.001; %Smalest depth value of the generation curve,eccept zero
                            stop=cell2mat(d_layer(Substrate_pos));
                            DP=90; % Nr of points for the generation profile (another three points will always be added to this value)
                            f=nthroot(stop/2/start,DP/2); % Calculation of a factor making the depth value increaing exponentially (derived from formula: d_layer/2=start*f^i where i=[0:1:steps/2])
                            depth=[0 start];
                            if lambda_batch<=960;
                                for ii=1:DP
                                    if ii<=DP/2
                                        depth(2+ii)=depth(2+ii-1)*f;
                                    else
                                        depth(2+ii)=depth(2+ii-1)+(depth(DP-ii+3)-depth(DP-ii+2)); % The point distribution calculated for the first half of the points is mirrord to the second half
                                    end
                                end
                                depth(length(depth)+1)=stop; % Final depth value is added to array so that the total array length is DP+3
                            else
                                depth=[0 2.^([0:DP+1].*log2(stop)/(DP+1))]; % Creates
                            end
                            
                            depth_eff=depth+sum(cell2mat(d_layer(1:Substrate_pos-1)));
                                                        
                            Gen_eff = histcounts(-real([A_tot_struct([A_tot_struct.layer_position]==Substrate_pos).abs_loc_eff]),depth_eff);
                            Gen_eff=[0 Gen_eff]'; % The Gen variable has always one value less then the depth variable. And the first value for PC1D needs to be 0/0.
                            Gen_eff=[depth'-depth(1) Gen_eff];
                                                        
                            % **** Writing generation profile in ASCII-format with .gen-suffix ****
                            if gen_create==1
                                h_Js=6.62607e-34;
                                c=299792458;
                                
                                ph_num=100*(lambda_batch/1e9)/(h_Js*c)/1e4; % specific nr of photons per cm² and J at a certain wavelength
                                cd(Gen_folder);
                                fid = fopen(['Gen_' num2str(lambda_batch) '.gen'], 'wt');
                                for i=1:length(Gen_eff(:,1))
                                    fprintf(fid, '%f\t' , [Gen_eff(i,1) sum(Gen_eff(1:i,2))*ph_num/sum(Gen_eff(1:end,2))]); % Dividing through the total sum is done to normalize the generation curve to 1
                                    fprintf(fid, '\n');
                                    Gen_PC1D(i,1:2)=[Gen_eff(i,1) sum(Gen_eff(1:i,2))*ph_num/sum(Gen_eff(1:end,2))];
                                end
                                fclose(fid);
                            end
                            assignin('base', ['Gen_PC1D_' num2str(lambda_batch)], Gen_PC1D)
                        end
                        
                        %% TESTING if sum of all processes equals 1
                        if abs(sum(R_esc(:,2)+A_fi(:,2)+A_bb(:,2)+A_AR(:,2)+A_1(:,2)+A_2(:,2)+A_3(:,2)+A_4(:,2)+T_tot(:,2)+error_tot(:,2)+interrupt_tot(:,2))/row-1)>0.0001
                            errordlg('SOME RAYS WERE LOST!')
                            R_esc(:,2)+A_fi(:,2)+A_bb(:,2)+A_AR(:,2)+A_1(:,2)+A_2(:,2)+A_3(:,2)+A_4(:,2)+T_tot(:,2)+error_tot(:,2)+interrupt_tot(:,2)
                            return
                        end                     
                        
                        %% Saving all variable and cleaning up workspace                        
                        save([datestr(now,formatOut) '_Results_batch_' num2str(batch) '_angle_' num2str(theta) '.mat']);
                        clearvars -except save_file_loc lambda T_all_batch batch alpha theta pyr_angle_glass pyr_angle_si d_glass AR
                                                
                    end
                end
            end
        end
    end
end
clear all


