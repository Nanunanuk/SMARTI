% nicht bei pwd mit der auswahl starten
% Link zu PC1D nur einmal eingeben
% EQE berechnen und plotten
% Beschreibung des codes

clear all
clc
[gen_file_loc] = uigetdir2(pwd,'Select folder with generation profiles');


q_C=1.6021766e-19;


for batch_cnt=1:length(gen_file_loc)
    cd('C:\Users\Nils\Documents\Nils Reiners\Simulation und Software\PC1D\PC1Dmod and cmd-PC1D v6.2.1');
    count=0;
    for lambda=[300:30:1200]
        count=count+1;
        wl(count)=lambda;
        %%%%%%%%%%%%%%%%%%  Konvertiere prm-Datei in txt-Datei   %%%%%%%%%%%%%%%%%%%%%%%%%%
        convert_str=['convert_prm_to_ascii.exe Pvcell_GenBatch_fixedLS.prm Pvcell_GenBatch_fixedLS.txt'];
        system(convert_str);
        fileformat='%s';
        
        %%%%%%%%%%%%%%%%%%  Gen_file einlesen   %%%%%%%%%%%%%%%%%%%%%%%%%%
        fileformat='%f%f';
        fileID=fopen([gen_file_loc{batch_cnt} '\Gen_' num2str(lambda) '.gen']);
        Daten1=textscan(fileID,fileformat, 'Delimiter', '');
        fclose(fileID);
        
        %%%%%%%%%%%%%%%%%%  Parameter txt-Datei einlesen   %%%%%%%%%%%%%%%%%%%%%%%%%%
        fileformat='%s';
        
        fileID=fopen('C:\Users\Nils\Documents\Nils Reiners\Simulation und Software\PC1D\PC1Dmod and cmd-PC1D v6.2.1\Pvcell_GenBatch_fixedLS.txt');
        Daten2=textscan(fileID,fileformat, 'Delimiter', '');
        fclose(fileID);
        
        %%%%%%%%%%%%%%%%%%  Parameter ändern   %%%%%%%%%%%%%%%%%%%%%%%%%%
        old_str_1=char(Daten2{1}(300));
        old_str_2=char(Daten2{1}(307));
        
        new_str=[gen_file_loc{batch_cnt} '\Gen_' num2str(lambda) '.gen'];
        
        Daten2{1}(300)=cellstr([old_str_1(1:20) new_str]);
        Daten2{1}(307)=cellstr([old_str_2(1:29) new_str]);
        
        %%%%%%%%%%%%%%%%%%  Neue Parameter txt-Datei schreiben   %%%%%%%%%%%%%%%%%%%%%%%%%%
        !del Pvcell_GenBatch_fixedLS.txt
        
        fid = fopen('Pvcell_GenBatch_fixedLS.txt', 'wt');
        for i=1:length(Daten2{1})
            fprintf(fid, '%s\n' , char(Daten2{1}(i)));
        end
        fclose(fid);
        
        %%%%%%%%%%%%%%%%%%  Konvertiere txt-Datei in prm-Datei   %%%%%%%%%%%%%%%%%%%%%%%%%%
        !del Pvcell_GenBatch_fixedLS.prm
        
        convert_str=['convert_ascii_to_prm.exe Pvcell_GenBatch_fixedLS.txt Pvcell_GenBatch_fixedLS.prm'];
        system(convert_str);
        
        %%%%%%%%%%%%%%%%%%  Simulation starten   %%%%%%%%%%%%%%%%%%%%%%%%%%
        [status,data_string]=system('cmd-pc1d6-2.exe Pvcell_GenBatch_fixedLS.prm')
        
        iv=str2num(data_string(38:end));
        curr(count)=iv(2);
        max_vals(count)=max(Daten1{2});
        
    end
    IQE=[wl' ((-curr/100/q_C)./max_vals)'];

    assignin('base',['IQE_' num2str(batch_cnt)], IQE);
    cd(gen_file_loc{batch_cnt});
    
%     load('A_Si.mat')
%     j=sum(q_C*IQE.*A_Si;);
%     save('j.mat','j');
    save('IQE.mat','IQE');
    plot(IQE(:,1),IQE(:,2)) % Strom wird transformiert in [mA/cm²] und max_vals hat die EInheit [1/cm²]
    hold on
    clear IQE
end





