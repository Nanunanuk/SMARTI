%|************************************************************************
%|SMARTI - Layer function
%|by Nils Reiners
%|This source code is free: you can redistribute it and/or modify it under the terms of the GNU General Public License (GPL) as published by the Free Software Foundation., version 3 of the License. Any redistribution of		 |% %|modified code must also be free under the GNU General Public License.																																							 |%
%|This source code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for |% %|more details. You should have received a copy of the GNU General Public License along with this source code. If not, see <http://www.gnu.org/licenses/>.																		 |%
%|************************************************************************

% 

function [T_top, T_bot, A,error_struct, interrupt]=layer(batch,theta_batch,j,w,layer_data,layer_position,grid, S, I_top, I_bot, plotting,xth_ray)

S=[S.S_top S.S_side S.S_bot];

%% Creating structs that are needed later
T_top=struct('nr',[],'wavelength',[],'start_loc',[],'dir_vector',[],'ray_col',[]);
T_bot=struct('nr',[],'wavelength',[],'start_loc',[],'dir_vector',[],'ray_col',[]);
A=struct('nr',[],'wavelength',[],'layer_position',[],'abs_loc',[],'abs_loc_eff',[],'ray_col',[]);
error_struct=struct('nr',[],'wavelength',[],'start_loc',[],'dir_vector',[],'ray_col',[],'layer',[],'nr_of_bounces',[]);
interrupt=struct('nr',[],'wavelength',[],'start_loc',[],'dir_vector',[],'ray_col',[],'layer',[]);

%% Merging all rays in one variable
if isempty([I_bot.nr])==1&&isempty([I_top.nr])==1
    return
elseif isempty([I_bot.nr])==1&&isempty([I_top.nr])~=1
    I=I_top;
elseif isempty([I_bot.nr])~=1&&isempty([I_top.nr])==1
    I=I_bot;
elseif isempty([I_bot.nr])~=1&&isempty([I_top.nr])~=1
    I=[I_top I_bot];
end

%% Initializing counting variables
abs_count=0;
trans_count_top=0;
trans_count_bot=0;
count_for=0;                                                                % This variable is needed for the process display
count_hundred=0;                                                            % This variable is needed for the process display

%% If the layer is set to 'Reflector' with reflectivity between 0 and 1 all rays are counted as absorbed
if length([layer_data(layer_position).layer_mat])==1 
    A=I;
    return
else
    %% Loop for every ray
    for l=1:length(I); % Loop for every ray
        count_for=count_for+1;
        
        %% STATUS DISPLAY
        if count_hundred<fix(count_for/100)*100
            count_hundred=fix(count_for/100)*100;
            Status=['Batch: ' num2str(batch) '; Theta: ' num2str(theta_batch) '; Lambda: ' num2str(I(l).wavelength) '; Repeat: ' num2str(j) '; Layer: ' num2str(layer_position) '; Ray_nr: ' num2str(count_hundred)]
        end
        
        
        %% Calculation of absorption coefficient alpha
        c=299792458; %Speed of light in vacuum [m/s]
        wavelength=I(l).wavelength; %nm
        f=c/(wavelength/10^9); %1/s
        omega=2*pi*f;
        mat_row=find(layer_data(layer_position).layer_mat==wavelength,1);
        k=layer_data(layer_position).layer_mat(mat_row,3);
        alpha=2*k*omega/c/100; % 1/cm
        
        %% Initializing aditional variables
        t=0;
        side_hit=0;
        prev_hit='';                                                       % Initialisation of a variable to store the previous value of 'top','side' or 'bot'
        
        %% Determination of the depth at which the ray will be absorbed with the help of a random number
        if alpha==0;
            xG=inf;
        else
            xG=log(1-randi(1000000)/1000000)/-alpha*1e4; % Absorption depth [µm] (Beer-Lambert-Law)
        end
        
        %% Maximum number of transitions
        max_transitions=1000;
        
        %% Preallocating rst_struct (improving performance)
        xx(1:length(S))={0};
        rst_struct_init=struct('nr',xx,'r',xx,'s',xx,'t',xx,'first_hit',xx,'pos',xx);
                
        %% Taking the starting point and direction vector from the input data
        p=I(l).start_loc;
        P=I(l).dir_vector;
        
        %% Starting bouncing of ray through layer
        count_while=0;
        while count_while<max_transitions % Loop of intersections per ray
            
            %% Checking if maximum number of transitions is reached
            %  => if yes, storring as interrupt
            if count_while==max_transitions-1;
                I_interrupt=I(l);
                I_interrupt.layer=layer_position;
                if isempty(interrupt(1).nr)==1
                    interrupt=I_interrupt;
                else
                    interrupt=[interrupt I_interrupt];
                end
                break
            end
            
            %% Counting number of intersections
            count_while=count_while+1;
            
            %% Variable for grid algorythm
            grid_refl=0;
            
            %% ***DETERMINING INTERSECTION***            
            [rst_struct hit]=intersect(p,P,S,rst_struct_init);
            
            %% Plotting if plotting is activated
            if plotting==1&&round(l/xth_ray)==l/xth_ray
                ray_col=I(l).ray_col;
                figure(1)
                hold on
                plot3([p(1) p(1)+rst_struct(hit).t*P(1)],[p(2) p(2)+rst_struct(hit).t*P(2)],[p(3) p(3)+rst_struct(hit).t*P(3)],'color',ray_col,'linewidth',2)
                hold on
            end
            
            %% If no hit could be found the cooresponding ray is stored in error_struct
            if hit==0
                I_error=I(l);
                I_error.layer=layer_position;                               % The information in which layer the error accured is added
                I_error.nr_of_bounces=count_while;                          % The information of how many bounces were counted before the error accured is added
                if isempty(error_struct(1).nr)==1
                    error_struct=I_error;
                else
                    error_struct=[error_struct I_error];
                end
                break
            end
            
            %% Total length of ray
            t=t+rst_struct(hit).t;
            
            %% Testing if the ray is hitting top, side or bottom multiple times successively
            if isequal(prev_hit,S(hit).pos)==1
                multi_hit=multi_hit+1;
            else
                multi_hit=0;
            end
            prev_hit=S(hit).pos;
            
            %% TOP, SIDE OR BOTTOM
            
            switch (S(hit).pos)
                %% Top of the unit cell is intersected
                case{'top'}
                    
                    %% Testing if ray is absorbed on way to intersection
                    if t>xG; % Absobed on path?
                        abs_count=abs_count+1;
                        A(abs_count).nr=I(l).nr;
                        A(abs_count).wavelength=I(l).wavelength;
                        A(abs_count).layer_position=layer_position;
                        p2=p+rst_struct(hit).t*P; % Calculating point of intersection
                        P12=p2-p; % Vector of the last transmission step
                        c=t-xG; % Rest length from absorption point so intersection
                        Pg=p+(norm(P12)-c)*P12/norm(P12); % Vector indicating absorption point
                        A(abs_count).abs_loc=Pg(3); % z component of absorption point
                        % A(abs_count).abs_loc_eff=Pg(3);
                        
                        
                        %% Calculating effective generation path (for PC1D simulation)
                        if Pg(3)>-sum([layer_data(1:layer_position-1).d_layer])-1/6*[layer_data(layer_position-1).p_bot]
                            if Pg(1)>Pg(2)
                                if Pg(2)>-1*Pg(1)+w
                                    ds=w-Pg(1);
                                else
                                    ds=Pg(2);
                                end
                            else
                                if Pg(2)>-1*Pg(1)+w
                                    ds=w-Pg(2);
                                else
                                    ds=Pg(1);
                                end
                            end
                            a=atand(layer_data(layer_position-1).p_bot/(0.5*w));
                            dz2=-(Pg(3)+sum([layer_data(1:layer_position-1).d_layer]));
                            if tand(ds/dz2)<a
                                dz1=Pg(3)+sum([layer_data(1:layer_position-1).d_layer]);
                                dx=ds-dz1/tand(a);
                                zeta=sind(a)*dx;
                                zeta=-sum([layer_data(1:layer_position-1).d_layer])-zeta;
                            else
                                zeta=sqrt(dz2^2+ds^2);
                                zeta=-sum([layer_data(1:layer_position-1).d_layer])-zeta;
                            end
                            A(abs_count).abs_loc_eff=zeta;
                        else
                            % A(abs_count).abs_loc=Pg(3); % If the absorption does not take place in the designated region beneath the pyramide
                            A(abs_count).abs_loc_eff=Pg(3);
                        end
                        
                        break
                    else
                        
                        %******* GRID **************
                        if grid.pos{layer_position}==1&&multi_hit==0 % testing if grid is implemented and if the texture is hit the first time
                            
                            rand_num2=randi(1000000)/1000000;
                            rand_num3=randi(1000000)/1000000;
                            % n1:
                            komplex_n1=layer_data(layer_position).layer_mat(mat_row,2)-layer_data(layer_position).layer_mat(mat_row,3)*i;
                            
                            if grid.share{1}>rand_num2 % Absorption in finger?
                                % n2:
                                komplex_n2=grid.mat{1}(mat_row,2)-grid.mat{1}(mat_row,3)*i;
                                % Reflection or absorption?
                                theta_rad=acos(dot(S(hit).N,-P));
                                [R,T_AR,A_AR]=AR_RAT(komplex_n1,komplex_n1,komplex_n2,wavelength,0,theta_rad);
                                if rand_num3<R % Ray is reflected! 
                                    p=p+rst_struct(hit).t*P;
                                    P=P+2*dot(-P,S(hit).N)*S(hit).N; % Angle of incidence equals angle of reflection
                                    grid_refl=1;
                                else % Ray is absorbed
                                    abs_count=abs_count+1;
                                    A(abs_count).nr=I(l).nr;
                                    A(abs_count).wavelength=I(l).wavelength;
                                    A(abs_count).layer_position=-1;
                                    break
                                end
                                
                            elseif grid.share{1}+grid.share{2}>rand_num2 % Absorption in busbar?
                                % n2:
                                komplex_n2=grid.mat{2}(mat_row,2)-grid.mat{2}(mat_row,3)*i;
                                % Reflection or absorption?
                                theta_rad=acos(dot(S(hit).N,-P));
                                [R,T_AR,A_AR]=AR_RAT(komplex_n1,komplex_n1,komplex_n2,wavelength,0,theta_rad);
                                if rand_num3<R % Ray is reflected! 
                                    p=p+rst_struct(hit).t*P;
                                    P=P+2*dot(-P,S(hit).N)*S(hit).N; % Angle of incidence equals angle of reflection
                                    grid_refl=1;
                                else % Ray is absorbed
                                    abs_count=abs_count+1;
                                    A(abs_count).nr=I(l).nr;
                                    A(abs_count).wavelength=I(l).wavelength;
                                    A(abs_count).layer_position=-2;
                                    break
                                end
                                
                            end
                        end
                        
                        %******* R,A,T material transition **************
                        
                        if grid_refl~=1 % R,A,T material transition is just calculated, if there has not been a reflection on the grid before
                            
                            % n1:
                            komplex_n1=layer_data(layer_position).layer_mat(mat_row,2)-layer_data(layer_position).layer_mat(mat_row,3)*i;
                            
                            % n_AR:
                            if length(layer_data(layer_position).AR_mat)<10
                                komplex_nAR=komplex_n1;
                                d_AR=0;
                            else
                                komplex_nAR=layer_data(layer_position).AR_mat(mat_row,2)-layer_data(layer_position).AR_mat(mat_row,3)*i;
                                d_AR=layer_data(layer_position).d_AR;
                            end
                            
                            % n2:
                            if layer_position==1
                                komplex_n2=1;
                            else
                                komplex_n2=layer_data(layer_position-1).layer_mat(mat_row,2)-layer_data(layer_position-1).layer_mat(mat_row,3)*i;
                            end
                            
                            % Reflection
                            theta_rad=acos(dot(S(hit).N,-P));
                            [R,T_AR,A_AR]=AR_RAT(komplex_n1,komplex_nAR,komplex_n2,wavelength,d_AR,theta_rad);
                            
                            rand_AR=randi(1000000)/1000000;
                            
                            if rand_AR<R % Reflection!
                                p=p+rst_struct(hit).t*P;
                                P=P+2*dot(-P,S(hit).N)*S(hit).N; % Angle of incidence equals angle of reflection
                            elseif rand_AR>R+T_AR % Absorption
                                abs_count=abs_count+1;
                                A(abs_count).nr=I(l).nr;
                                A(abs_count).wavelength=I(l).wavelength;
                                A(abs_count).layer_position=0;
                                break
                            else % Transmission
                                trans_count_top=trans_count_top+1;
                                T_top(trans_count_top).nr=I(l).nr;
                                T_top(trans_count_top).wavelength=I(l).wavelength;
                                T_top(trans_count_top).start_loc=p+rst_struct(hit).t*P;
                                T_top(trans_count_top).dir_vector=snell(layer_data,layer_position,layer_position-1,wavelength,P,S(hit).N);
                                T_top(trans_count_top).ray_col=I(l).ray_col;
                                break
                            end
                            
                        end
                    end
                    side_hit=0; % Resetting the side_hit counter for the "Booster"
                    %% Side of the unit cell is intersected
                case{'side'}
                    % Absorption?
                    if t>xG;
                        abs_count=abs_count+1;
                        A(abs_count).nr=I(l).nr;
                        A(abs_count).wavelength=I(l).wavelength;
                        A(abs_count).layer_position=layer_position;
                        p2=p+rst_struct(hit).t*P;
                        P12=p2-p;
                        c=t-xG;
                        Pg=p+(norm(P12)-c)*P12/norm(P12);
                        A(abs_count).abs_loc=Pg(3);
                        % A(abs_count).abs_loc_eff=Pg(3);
                        
                        if Pg(3)>-sum([layer_data(1:layer_position-1).d_layer])-1/6*[layer_data(layer_position-1).p_bot]
                            if Pg(1)>Pg(2)
                                if Pg(2)>-1*Pg(1)+w
                                    ds=w-Pg(1);
                                else
                                    ds=Pg(2);
                                end
                            else
                                if Pg(2)>-1*Pg(1)+w
                                    ds=w-Pg(2);
                                else
                                    ds=Pg(1);
                                end
                            end
                            a=atand(layer_data(layer_position-1).p_bot/(0.5*w));
                            dz2=-(Pg(3)+sum([layer_data(1:layer_position-1).d_layer]));
                            if tand(ds/dz2)<a
                                dz1=Pg(3)+sum([layer_data(1:layer_position-1).d_layer]);
                                dx=ds-dz1/tand(a);
                                zeta=sind(a)*dx;
                                zeta=-sum([layer_data(1:layer_position-1).d_layer])-zeta;
                            else
                                zeta=sqrt(dz2^2+ds^2);
                                zeta=-sum([layer_data(1:layer_position-1).d_layer])-zeta;
                            end
                            A(abs_count).abs_loc_eff=zeta;
                        else
                            %A(abs_count).abs_loc=Pg(3);
                            A(abs_count).abs_loc_eff=Pg(3);
                        end
                        break
                    else
                        
                        
                        
                        %% SPEEDING UP OF THE TRANSITION THROUGH A THICK LAYER
                        side_hit=side_hit+1;
                        Booster=1; %Activate algorithm?
                        if side_hit==2&&Booster==1
                                                                                       
                                if P(3)<=0 % checking if the ray is pointing down or up
                                    screenpos=-sum([layer_data(1:layer_position).d_layer])+layer_data(layer_position).p_bot+0.1; % Calculating the z componant of a horizontal layer above the bottom texture
                                    S_botscreen=struct('p',[w/2;w/2;screenpos],'v1',[1;0;0],'v2',[0;1;0],'N',[0;0;1],'pos','bot_scr'); % Determining the representation of the layer for the intersect algorythm
                                    rst_struct_init=struct('nr',0,'r',0,'s',0,'t',0,'first_hit',0,'pos',0); % Resetting rst_struct
                                    [rst2]=intersect(p,P,S_botscreen,rst_struct_init); % Calculating rst for the ray and the layer
                                    h0=abs(screenpos)-abs(p(3));
                                elseif P(3)>0 % Upwards
                                    if layer_position==1
                                        screenpos=-0.1;
                                    else
                                        screenpos=-sum([layer_data(1:layer_position-1).d_layer])-0.1;
                                    end
                                    S_topscreen=struct('p',[w/2;w/2;screenpos],'v1',[0;1;0],'v2',[1;0;0],'N',[0;0;-1],'pos','top_scr');
                                    rst_struct_init=struct('nr',0,'r',0,'s',0,'t',0,'first_hit',0,'pos',0);
                                    [rst2]=intersect(p,P,S_topscreen,rst_struct_init);
                                    h0=abs(p(3))-abs(screenpos);
                                end
                                t=t-rst_struct(hit).t+rst2.t;
                                if t>xG;
                                    abs_count=abs_count+1;
                                    A(abs_count).nr=I(l).nr;
                                    A(abs_count).wavelength=I(l).wavelength;
                                    A(abs_count).layer_position=layer_position;
                                    p2=p+rst2.t*P;
                                    P12=p2-p; %Distance between hitpoint on screen and start position at side
                                    c=t-xG;
                                    Pg=p+(norm(P12)-c)*P12/norm(P12);
                                    A(abs_count).abs_loc=Pg(3);
                                    A(abs_count).abs_loc_eff=Pg(3);
                                    break
                                else
                                    r=sqrt(rst2.t^2-h0^2);
                                    pneu=p+rst2.t*P;
                                    if pneu(1)>=0&&pneu(2)>=0
                                        p=[abs(pneu(1));abs(pneu(2));screenpos]-[floor(abs(pneu(1))/w)*w;floor(abs(pneu(2))/w)*w;0];
                                    elseif pneu(1)<0&&pneu(2)>=0
                                        p=[w;0;0]-([abs(pneu(1));abs(pneu(2));screenpos]-[floor(abs(pneu(1))/w)*w;floor(abs(pneu(2))/w)*w;0]);
                                        p=[abs(p(1));abs(p(2));-p(3)];
                                    elseif pneu(1)<0&&pneu(2)<0
                                        p=[w;w;0]-([abs(pneu(1));abs(pneu(2));screenpos]-[floor(abs(pneu(1))/w)*w;floor(abs(pneu(2))/w)*w;0]);
                                        p=[abs(p(1));abs(p(2));-p(3)];
                                    elseif pneu(1)>=0&&pneu(2)<0
                                        p=[0;w;0]-([abs(pneu(1));abs(pneu(2));screenpos]-[floor(abs(pneu(1))/w)*w;floor(abs(pneu(2))/w)*w;0]);
                                        p=[abs(p(1));abs(p(2));-p(3)];
                                    end                                    
                                end
                                                                                                                   
                        else % If the side is hit for the first time the booster function is not active
                            h=p+rst_struct(hit).t*P; %Here the switch from one side of the unit cell to the other is calculated
                            if h(1)<0.000001&&h(1)>-0.000001
                                p=p+rst_struct(hit).t*P+[w;0;0];
                            elseif h(1)<w+0.000001&&h(1)>w-0.000001
                                p=p+rst_struct(hit).t*P-[w;0;0];
                            elseif h(2)<0.000001&&h(2)>-0.000001
                                p=p+rst_struct(hit).t*P+[0;w;0];
                            else
                                p=p+rst_struct(hit).t*P-[0;w;0];
                            end
                        end
                    end
                    %% Bottom of the unit cell is intersected
                case{'bot'}
                    %% Absorption
                    if t>xG;
                        abs_count=abs_count+1;
                        A(abs_count).nr=I(l).nr;
                        A(abs_count).wavelength=I(l).wavelength;
                        A(abs_count).layer_position=layer_position;
                        p2=p+rst_struct(hit).t*P;
                        P12=p2-p;
                        c=t-xG;
                        Pg=p+(norm(P12)-c)*P12/norm(P12);
                        A(abs_count).abs_loc=Pg(3);
                        %A(abs_count).abs_loc_eff=Pg(3);
                        
                        if Pg(3)>-sum([layer_data(1:layer_position-1).d_layer])-1/6*[layer_data(layer_position-1).p_bot]
                            if Pg(1)>Pg(2)
                                if Pg(2)>-1*Pg(1)+w
                                    ds=w-Pg(1);
                                else
                                    ds=Pg(2);
                                end
                            else
                                if Pg(2)>-1*Pg(1)+w
                                    ds=w-Pg(2);
                                else
                                    ds=Pg(1);
                                end
                            end
                            a=atand(layer_data(layer_position-1).p_bot/(0.5*w));
                            dz2=-(Pg(3)+sum([layer_data(1:layer_position-1).d_layer]));
                            if tand(ds/dz2)<a
                                dz1=Pg(3)+sum([layer_data(1:layer_position-1).d_layer]);
                                dx=ds-dz1/tand(a);
                                zeta=sind(a)*dx;
                                zeta=-sum([layer_data(1:layer_position-1).d_layer])-zeta;
                            else
                                zeta=sqrt(dz2^2+ds^2);
                                zeta=-sum([layer_data(1:layer_position-1).d_layer])-zeta;
                            end
                            A(abs_count).abs_loc_eff=zeta;
                        else
                            %A(abs_count).abs_loc=Pg(3);
                            A(abs_count).abs_loc_eff=Pg(3);
                        end
                        break
                        
                    else
                        
                        %******* GRID **************
                        if layer_position<size(layer_data,2) % Testing if ray is in the last layer
                            if grid.pos{layer_position+1}==1&&multi_hit==0 % testing if grid is implemented and if the texture is hit the first time
                                
                                rand_num2=randi(1000000)/1000000;
                                rand_num3=randi(1000000)/1000000;
                                % n1:
                                komplex_n1=layer_data(layer_position).layer_mat(mat_row,2)-layer_data(layer_position).layer_mat(mat_row,3)*i;
                                
                                if rand_num2<grid.share{1} % Absorption in finger?
                                    % n2:
                                    komplex_n2=grid.mat{1}(mat_row,2)-grid.mat{1}(mat_row,3)*i;
                                    % Reflection
                                    theta_rad=acos(dot(S(hit).N,-P));
                                    [R,T_AR,A_AR]=AR_RAT(komplex_n1,komplex_n1,komplex_n2,wavelength,0,theta_rad);
                                    if rand_num3<R
                                        p=p+rst_struct(hit).t*P;
                                        p(3)=-sum([layer_data(1:layer_position).d_layer])+layer_data(layer_position).p_bot;
                                        P=P+2*dot(-P,[0;0;1])*[0;0;1];
                                        grid_refl=1;
                                    else
                                        abs_count=abs_count+1;
                                        A(abs_count).nr=I(l).nr;
                                        A(abs_count).wavelength=I(l).wavelength;
                                        A(abs_count).layer_position=-1;
                                        break
                                    end
                                    
                                elseif rand_num2<grid.share{1}+grid.share{2} % Absorption in busbar?
                                    % n2:
                                    komplex_n2=grid.mat{2}(mat_row,2)-grid.mat{2}(mat_row,3)*i;
                                    % Reflection
                                    theta_rad=acos(dot(S(hit).N,-P));
                                    [R,T_AR,A_AR]=AR_RAT(komplex_n1,komplex_n1,komplex_n2,wavelength,0,theta_rad);
                                    if rand_num3<R
                                        p=p+rst_struct(hit).t*P;
                                        p(3)=-sum([layer_data(1:layer_position).d_layer])+layer_data(layer_position).p_bot;
                                        P=P+2*dot(-P,[0;0;1])*[0;0;1];
                                        grid_refl=1;
                                    else
                                        abs_count=abs_count+1;
                                        A(abs_count).nr=I(l).nr;
                                        A(abs_count).wavelength=I(l).wavelength;
                                        A(abs_count).layer_position=-2;
                                        break
                                    end
                                    
                                end
                            end
                        
                        
                        %******* R,A,T material transition **************
                        
                        if grid_refl~=1 % R,A,T material transition is just calculated, if there has not been a reflection on the grid before
                            
                            if length([layer_data(layer_position+1).layer_mat])==1  % Testing if reflector is set for final layer
                                R=[layer_data(layer_position+1).layer_mat];         % The reflectivity is set to the value that is written in "layer_mat"
                                T_AR=1;                                             % The transmission value of the AR-coating is set to 1 as no coating is existing
                            else
                                % n1:
                                komplex_n1=layer_data(layer_position).layer_mat(mat_row,2)-layer_data(layer_position).layer_mat(mat_row,3)*i;
                                
                                % n_AR:
                                if layer_position+1>length(layer_data)
                                    komplex_nAR=komplex_n1;
                                    d_AR=0;
                                elseif length(layer_data(layer_position+1).AR_mat)<10
                                    komplex_nAR=komplex_n1;
                                    d_AR=0;
                                else
                                    komplex_nAR=layer_data(layer_position+1).AR_mat(mat_row,2)-layer_data(layer_position+1).AR_mat(mat_row,3)*i;
                                    d_AR=layer_data(layer_position+1).d_AR;
                                end
                                
                                % n2:
                                if layer_position==4
                                    komplex_n2=1;
                                else
                                    komplex_n2=layer_data(layer_position+1).layer_mat(mat_row,2)-layer_data(layer_position+1).layer_mat(mat_row,3)*i;
                                end
                                
                                % Reflection
                                theta_rad=acos(dot(S(hit).N,-P));
                                [R,T_AR,A_AR]=AR_RAT(komplex_n1,komplex_nAR,komplex_n2,wavelength,d_AR,theta_rad);
                            end
                            
                            rand_AR=randi(1000000)/1000000;
                            if rand_AR<=R % Reflected
                                p=p+rst_struct(hit).t*P;
                                rand_lambert=randi(1000000)/1000000;
                                if layer_data(layer_position).lambert<rand_lambert
                                    P=P+2*dot(-P,S(hit).N)*S(hit).N; % Angle of incidence equals angle of reflection
                                else
                                    P=lambertsch(S(hit).N);                 % Gives back any diffuse reflected ray on surface S(hit)
                                end
                            elseif rand_AR>R+T_AR % Absorbed
                                abs_count=abs_count+1;
                                A(abs_count).nr=I(l).nr;
                                A(abs_count).wavelength=I(l).wavelength;
                                A(abs_count).layer_position=0;
                                break
                            else
                                % Transmitted
                                trans_count_bot=trans_count_bot+1;
                                T_bot(trans_count_bot).nr=I(l).nr;
                                T_bot(trans_count_bot).wavelength=I(l).wavelength;
                                T_bot(trans_count_bot).start_loc=p+rst_struct(hit).t*P;
                                if length([layer_data(layer_position+1).layer_mat])==1 % Testing if reflector is set for final layer
                                    T_bot(trans_count_bot).dir_vector=I(l).dir_vector;
                                else
                                    T_bot(trans_count_bot).dir_vector=snell(layer_data,layer_position,layer_position+1,wavelength,P,S(hit).N);
                                end
                                T_bot(trans_count_bot).ray_col=I(l).ray_col;
                                break
                            end
                        end
                        else % If the ray is in the last layer the rays are all counted as transmitted if they hit the backside
                                % Transmitted
                                trans_count_bot=trans_count_bot+1;
                                T_bot(trans_count_bot).nr=I(l).nr;
                                T_bot(trans_count_bot).wavelength=I(l).wavelength;
                                T_bot(trans_count_bot).start_loc=p+rst_struct(hit).t*P; 
                                T_bot(trans_count_bot).dir_vector=I(l).dir_vector;
                                T_bot(trans_count_bot).ray_col=I(l).ray_col;
                                break
                        end
                        end
                    side_hit=0; % Resetting the side_hit counter for the "Booster"
            end
        end
    end
end
