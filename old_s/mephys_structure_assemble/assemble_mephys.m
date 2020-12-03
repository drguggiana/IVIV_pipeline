%% L23 Morpho
filename=uipickfiles('FilterSpec','D:\Postdoc_Margrie\Projects\L23\Morph_ephys')%pathname, you need uipickfiles function
load(char(filename));%load mat file

data_morpho=com_TMDidx;
data_morpho(:,[9 12 21 24 26 27 28])=[];
pia_morpho=com_TMDidx(:,26);
cellID_morpho=com_TMDidx(:,28);
%Cluster analysis performed via TMD in Phython
idx_morpho=com_TMDidx(:,27);
%% load traces of morpho
filename=uipickfiles('FilterSpec','D:\Postdoc_Margrie\Projects\L23\Morph_ephys')%pathname, you need uipickfiles function
load(char(filename));%load mat file

%% L23 EPHYS
filename=uipickfiles('FilterSpec','D:\Postdoc_Margrie\Projects\L23\Morph_ephys\')%pathname, you need uipickfiles function
load(char(filename));%load mat file

%% Save location

adata_dir='C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Full_structure\mephys_structure'

%% EPHYS Clustering 
%Prepare data 
data_ephys=reshape([ephys_thesis{:,2:19}],[137,18])
pia_ephys=reshape([ephys_thesis{:,20}],[137,1]);
cellID_ephys=reshape([ephys_thesis{:,21}],[137,1]);

%% Assemblelarge structure
[~,ia,ib] =intersect(cellID_ephys,cellID_morpho);

%% Morpho
temp=data_morpho;
temp(ib,:)=[];
kj=1:189
kj(ib)=[];
for i=1:length(kj)
  cnames{i}=[traces{kj(i),4}];
end
 
for i=1:137
enames{i}=ephys_thesis{i,1}
end

menames=[enames cnames];

%% 
for i=1:293
mephys(i).cell_name=[menames{i}];
end
%% 

for i=1:137

mephys(i).id=[ephys_thesis{i,21}];   
mephys(i).Vmin=[ephys_thesis{i,2}];  
mephys(i).Vpeak=[ephys_thesis{i,3}];  
mephys(i).Vthresh=[ephys_thesis{i,4}];  
mephys(i).Vslope=[ephys_thesis{i,5}];  
mephys(i).Vhalf=[ephys_thesis{i,6}];  
mephys(i).Vamp=[ephys_thesis{i,7}];  
mephys(i).AHP=[ephys_thesis{i,8}]; 
mephys(i).APrise=[ephys_thesis{i,9}]; 
mephys(i).APfall=[ephys_thesis{i,10}]; 
mephys(i).APbwidth=[ephys_thesis{i,11}]; 
mephys(i).APhwidth=[ephys_thesis{i,12}]; 
mephys(i).APlate=[ephys_thesis{i,13}]; 
mephys(i).Vrest=[ephys_thesis{i,14}]; 
mephys(i).tau=[ephys_thesis{i,15}]; 
mephys(i).Rin=[ephys_thesis{i,16}]; 
mephys(i).Sag=[ephys_thesis{i,17}]; 
mephys(i).Rheo=[ephys_thesis{i,18}]; 
mephys(i).APfreq=[ephys_thesis{i,19}];
mephys(i).pia=[ephys_thesis{i,20}];
end

%% 
temp2=com_TMDidx;
temp2(ib,:)=[];
lf=138:293;

for i=1:156
mephys(lf(i)).id=temp2(i,28); 
mephys(lf(i)).RDA=temp(i,1); 
mephys(lf(i)).LA=temp(i,2); 
mephys(lf(i)).PLA=temp(i,3); 
mephys(lf(i)).BPA=temp(i,4);
mephys(lf(i)).BOA=temp(i,5);
mephys(lf(i)).BLA=temp(i,6);
mephys(lf(i)).PLA=temp(i,7);
mephys(lf(i)).WHA=temp(i,8);
mephys(lf(i)).XSA=temp(i,9);
mephys(lf(i)).YSA=temp(i,10);

mephys(lf(i)).RDB=temp(i,11); 
mephys(lf(i)).LB=temp(i,12); 
mephys(lf(i)).PLB=temp(i,13); 
mephys(lf(i)).BPB=temp(i,14);
mephys(lf(i)).BOB=temp(i,15);
mephys(lf(i)).BLB=temp(i,16);
mephys(lf(i)).PLB=temp(i,17);
mephys(lf(i)).WHB=temp(i,18);
mephys(lf(i)).XSB=temp(i,19);
mephys(lf(i)).YSB=temp(i,20);
mephys(lf(i)).NB=temp(i,21);
mephys(lf(i)).TMD_idx=temp2(i,27);

mephys(lf(i)).pia=temp2(i,26);

mephys(lf(i)).tr_apical=traces{kj(i),1};
mephys(lf(i)).tr_basal=traces{kj(i),2};
mephys(lf(i)).tr_soma=traces{kj(i),3};
end
%% 
for i=1:length(ib)
    mephys(ia(i)).RDA=data_morpho(ib(i),1); 
    mephys(ia(i)).LA=data_morpho(ib(i),2); 
    mephys(ia(i)).PLA=data_morpho(ib(i),3); 
    mephys(ia(i)).BPA=data_morpho(ib(i),4); 
    mephys(ia(i)).BOA=data_morpho(ib(i),5); 
    mephys(ia(i)).BLA=data_morpho(ib(i),6); 
    mephys(ia(i)).PLA=data_morpho(ib(i),7); 
    mephys(ia(i)).WHA=data_morpho(ib(i),8); 
    mephys(ia(i)).XSA=data_morpho(ib(i),9); 
    mephys(ia(i)).YSA=data_morpho(ib(i),10); 
    
      mephys(ia(i)).RDB=data_morpho(ib(i),11); 
    mephys(ia(i)).LB=data_morpho(ib(i),12); 
    mephys(ia(i)).PLB=data_morpho(ib(i),13); 
    mephys(ia(i)).BPB=data_morpho(ib(i),14); 
    mephys(ia(i)).BOB=data_morpho(ib(i),15); 
    mephys(ia(i)).BLB=data_morpho(ib(i),16); 
    mephys(ia(i)).PLB=data_morpho(ib(i),17); 
    mephys(ia(i)).WHB=data_morpho(ib(i),18); 
    mephys(ia(i)).XSB=data_morpho(ib(i),19); 
    mephys(ia(i)).YSB=data_morpho(ib(i),20); 
    mephys(ia(i)).NB=data_morpho(ib(i),21); 
    
     mephys(ia(i)).TMD_idx=com_TMDidx(ib(i),27); 
     
    mephys(ia(i)).tr_apical=traces{ib(i),1};
mephys(ia(i)).tr_basal=traces{ib(i),2};
mephys(ia(i)).tr_soma=traces{ib(i),3};
end


%% 
% 
% 
% for i=1:length(ib)
% % mephys(i).id=[]; 
% % mephys(i).Vmin=[];  
% % mephys(i).Vpeak=[];  
% % mephys(i).Vthresh=[];  
% % mephys(i).Vslope=[];  
% % mephys(i).Vhalf=[];  
% % mephys(i).Vamp=[];  
% % mephys(i).AHP=[]; 
% % mephys(i).APrise=[]; 
% % mephys(i).APfall=[]; 
% % mephys(i).APbwidth=[]; 
% % mephys(i).APhwidth=[]; 
% % mephys(i).APlate=[]; 
% % mephys(i).Vrest=[]; 
% % mephys(i).tau=[]; 
% % mephys(i).Rin=[]; 
% % mephys(i).Sag=[]; 
% % mephys(i).Rheo=[]; 
% % mephys(i).APfreq=[];
% % 
% %     mephys(ib(i)).RDA=[]; 
% %     mephys(ib(i)).LA=[]; 
% %     mephys(ib(i)).PLA=[]; 
% %     mephys(ib(i)).BPA=[]; 
% %     mephys(ib(i)).BOA=[]; 
% %     mephys(ib(i)).BLA=[]; 
% %     mephys(ib(i)).PLA=[]; 
% %     mephys(ib(i)).WHA=[]; 
% %     mephys(ib(i)).XSA=[]; 
% %     mephys(ib(i)).YSA=[]; 
% %     
% %       mephys(ib(i)).RDB=[]; 
% %     mephys(ib(i)).LB=[]; 
% %     mephys(ib(i)).PLB=[]; 
% %     mephys(ib(i)).BPB=[]; 
% %     mephys(ib(i)).BOB=[]; 
% %     mephys(ib(i)).BLB=[]; 
% %     mephys(ib(i)).PLB=[]; 
% %     mephys(ib(i)).WHB=[]; 
% %     mephys(ib(i)).XSB=[]; 
% %     mephys(ib(i)).YSB=[]; 
% %     mephys(ib(i)).NB=[]; 
% %     
%     mephys(ib(i)).tr_apical=traces{ia(i),1};
% mephys(ib(i)).tr_basal=traces{ia(i),2};
% mephys(ib(i)).tr_soma=traces{ia(i),3};
% 
% mephys(i).pia=[];
% end
%% 
%% Save
  cd(adata_dir);
        FileName=['mephys_struct'];
        %         save(FileName,'-struct','LGN','-v7.3');
        save(FileName,'mephys','-v7.3');
        disp('FILE SAVED');
