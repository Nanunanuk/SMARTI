function load_data(Mat)

%% Loading radiation data and material parameters
curr_path=mfilename('fullpath');
[a b c]=xlsread([curr_path(1:length(curr_path)-9) 'Data\All_data_PV-Lighthouse.xlsx']);


for n = 1:size(Mat,1)
    row1=find(strcmp(Mat{n,1},b(1,:)));
    row2=find(strcmp(Mat{n,2},b(2,:)));
    row3=find(strcmp(Mat{n,3},b(3,:)));
    
    x=row1(ismember(row1,row2));
    data_col=x(ismember(x,row3));
    
    wl_n_k=[a(:,1) a(:,data_col) a(:,data_col+1)];
    assignin('base', [Mat{n}], wl_n_k);
end