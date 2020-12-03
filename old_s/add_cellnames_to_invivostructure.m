 
%% LOAD in vivo data structure L23PC data structure, you need uipickfiles function
out_dir='C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\structures_invivo'
directory=out_dir;% use cobined date structure named Data_SWMF_combined_qualitymoph atm:191804
filename=uipickfiles('FilterSpec',directory)%pathname, you need uipickfiles function
load(char(filename));%load mat file

%% Excel sheet directory
ExpXls            = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\L23PC_project_all_invivo.xlsx';%directory where excel batch file is located;change accordingly
%%   Read out excel info
batchopt          = parseExperimentsXls_L23(ExpXls);%calls the nested function parseExperimentsXls_dLGN and considers the user flag (1 or 0)
nummice           = length(batchopt.mouse);%length of experiments to be analyzed
%% 
%cell2=cell(70,1)
adder=1;%counting variable
for i=1:length(L23_PC)
   temp=[]
    temp=L23_PC(i).ivivROI;
   k=[];
    cellname=cell(length(temp),1)
  for k=1:length(temp)
   cellname{k,:}=[char(batchopt.exp_invitro{i}) num2str(temp(k),'%04.f')];
  end 
  
  tr(adder).cellname=cellname
adder=adder+1;
end

%% 
%% 
cellnames_all=cat(1,tr(:).cellname);
