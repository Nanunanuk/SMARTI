function [S]=geometry_cube_pyramide(pos,prev,w,layer_data, plotting, xth_ray)
    
%% Determining intersection points for geometry with previous layer
if pos>1;
    
    S_top=struct('p',[],'v1',[],'v2',[],'N',[],'pos',[]);
    S_top=prev.S_bot;
    
    for n=1:length(S_top)
        S_top(n).N=-S_top(n).N;
        S_top(n).pos='top';
    end
    
else
    
    %% Top
       
    p1=[0; 0; 0];
    p2=[w; 0; 0];
    p3=[w; w; 0];
    p4=[0; w; 0];
    
    M1=[p1 p2 p4];
    M2=[p3 p4 p2];
    
    S_top=struct('p',[],'v1',[],'v2',[],'N',[],'pos',[]);
    
    S_top(1).p =p1;          
    S_top(1).v1=p2-p1;       
    S_top(1).v2=p4-p1;       
    
    S_top(2).p =p3;          
    S_top(2).v1=p4-p3;       
    S_top(2).v2=p2-p3;       
    
    if plotting==1
    figure(1)
    hold on
    fill3(M1(1,:), M1(2,:), M1(3,:),'r')
    hold on
    alpha(0.3)
    fill3(M2(1,:), M2(2,:), M2(3,:),'r')
    alpha(0.3)
    
    xlabel('x')
    ylabel('y')
    end
    
    for n=1:2
        
        S_top(n).pos='top';
        S_top(n).N=-cross(S_top(n).v1,S_top(n).v2)/norm(cross(S_top(n).v1,S_top(n).v2));
    end
    
end


%% Side

S_side=struct('p',[],'v1',[],'v2',[],'N',[],'pos',[]);

    d_layer_prev=sum([layer_data(1:find([layer_data(:).layer_position]==pos-1)).d_layer]);
    d_layer=layer_data(find([layer_data(:).layer_position]==pos)).d_layer;

    p1=[0; 0; 0-d_layer_prev];
    p3=[w; 0; 0-d_layer_prev];
    p5=[w; w; 0-d_layer_prev];
    p7=[0; w; 0-d_layer_prev];
    
    p2=[w; 0; 0-d_layer-d_layer_prev];
    p4=[w; w; 0-d_layer-d_layer_prev];
    p6=[0; w; 0-d_layer-d_layer_prev];
    p8=[0; 0; 0-d_layer-d_layer_prev];
    
    M1=[p2 p1 p8];
    M2=[p2 p3 p1];
    M3=[p4 p3 p2];
    M4=[p4 p5 p3];
    
    M5=[p6 p5 p4];
    M6=[p6 p7 p5];
    M7=[p8 p7 p6];
    M8=[p8 p1 p7];
    
    if plotting==1
    figure(1)
    fill3(M1(1,:), M1(2,:), M1(3,:),'b')
    hold on
    alpha(0.1)
    fill3(M2(1,:), M2(2,:), M2(3,:),'b')
    alpha(0.1)
    fill3(M3(1,:), M3(2,:), M3(3,:),'b')
    alpha(0.1)
    fill3(M4(1,:), M4(2,:), M4(3,:),'b')
    alpha(0.1)
    
    fill3(M5(1,:), M5(2,:), M5(3,:),'b')
    alpha(0.1)
    fill3(M6(1,:), M6(2,:), M6(3,:),'b')
    alpha(0.1)
    fill3(M7(1,:), M7(2,:), M7(3,:),'b')
    alpha(0.1)
    fill3(M8(1,:), M8(2,:), M8(3,:),'b')
    alpha(0.3)
    end
    
    S_side(1).p =p2;          %2
    S_side(1).v1=p1-p2;       %1-2
    S_side(1).v2=p8-p2;       %8-2
    
    S_side(2).p =p2;          %2
    S_side(2).v1=p3-p2;       %3-2
    S_side(2).v2=p1-p2;       %1-2
    
    S_side(3).p =p4;          %4
    S_side(3).v1=p3-p4;       %3-4
    S_side(3).v2=p2-p4;       %2-4
    
    S_side(4).p =p4;          %4
    S_side(4).v1=p5-p4;       %5-4
    S_side(4).v2=p3-p4;       %3-4
    
    S_side(5).p =p6;          %6
    S_side(5).v1=p5-p6;       %5-6
    S_side(5).v2=p4-p6;       %4-6
    
    S_side(6).p =p6;          %6
    S_side(6).v1=p7-p6;       %7-6
    S_side(6).v2=p5-p6;       %5-6
    
    S_side(7).p =p8;          %8
    S_side(7).v1=p7-p8;       %7-8
    S_side(7).v2=p6-p8;       %6-8
    
    S_side(8).p =p8;          %8
    S_side(8).v1=p1-p8;       %1-8
    S_side(8).v2=p7-p8;       %7-8
    
    
for n=1:8
    S_side(n).N=-cross(S_side(n).v1,S_side(n).v2)/norm(cross(S_side(n).v1,S_side(n).v2));
    S_side(n).pos='side';
end

%% Bottom

    p_bot=layer_data(find([layer_data(:).layer_position]==pos)).p_bot;

    p1=[w/2; w/2; 0-d_layer-d_layer_prev+p_bot];
    p2=[0; 0; 0-d_layer-d_layer_prev];
    p3=[w; 0; 0-d_layer-d_layer_prev];
    p4=[w; w; 0-d_layer-d_layer_prev];
    p5=[0; w; 0-d_layer-d_layer_prev];
    
    M1=[p1 p2 p3];
    M2=[p1 p3 p4];
    M3=[p1 p4 p5];
    M4=[p1 p5 p2];
    
    S_bot=struct('p',[],'v1',[],'v2',[],'N',[],'pos',[]);
    
    S_bot(1).p =p1;          %1
    S_bot(1).v1=p2-p1;       %2-1
    S_bot(1).v2=p3-p1;       %3-1
    
    S_bot(2).p =p1;          %1
    S_bot(2).v1=p3-p1;       %3-1
    S_bot(2).v2=p4-p1;       %4-1
    
    S_bot(3).p =p1;          %1
    S_bot(3).v1=p4-p1;       %4-1
    S_bot(3).v2=p5-p1;       %5-1
    
    S_bot(4).p =p1;          %1
    S_bot(4).v1=p5-p1;       %5-1
    S_bot(4).v2=p2-p1;       %2-1
    
    if plotting==1
    figure(1)
    fill3(M1(1,:), M1(2,:), M1(3,:),'r')
    alpha(0.3)
    fill3(M2(1,:), M2(2,:), M2(3,:),'r')
    alpha(0.3)
    fill3(M3(1,:), M3(2,:), M3(3,:),'r')
    alpha(0.3)
    fill3(M4(1,:), M4(2,:), M4(3,:),'r')
    alpha(0.3)
    
    xlabel('x')
    ylabel('y')
    set(gca,'DataAspectRatio',[1 1 1])
    end
       
    for n=1:4
        
        S_bot(n).pos='bot';
        S_bot(n).N=cross(S_bot(n).v1,S_bot(n).v2)/norm(cross(S_bot(n).v1,S_bot(n).v2));
    end

S.S_top=S_top;
S.S_side=S_side;
S.S_bot=S_bot;