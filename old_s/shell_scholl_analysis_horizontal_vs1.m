function output=shell_scholl_analysis_horizontal_vs1(cyl_coord,sampling,extend,tree_type)
%VS 07/01/21
%input:
%cyl_coord - matrix with start and endpoint of each dendritic segment in a
%row from cyl function in trees toolbox
%sampling - sampling of x coordinate
%extend - interval to be covered
%tree_type= 1: apical, 2: basal
%output: see below

%example: shell_scholl_analysis_horizontal_vs1(cyl_coord,5,300)

figure;
hold
a=0
tempScholl=[];

for i=1:length(cyl_coord.apical)
    i
%     dataAll = cyl_coord.apical{1,i}; %for apical
    %dataAll = cyl_coord.apical{1,i}; %for apical
    if tree_type==1
    dataAll = cyl_coord.apical{1,i}; %for apical
    else tree_type==2
    dataAll = cyl_coord.basal{1,i}; %for basal
    end
    if ~isnan(dataAll)
        a=a+1;
        data=dataAll(:,1:2);
        output=scholl_analysis_horizontal_vs1(data,sampling,extend);
        tempScholl=[tempScholl; output.resampleScholl];
        plot(output.resampleX,output.resampleScholl,'blue')
    end
end
a
size(tempScholl)
Scholl=mean(tempScholl,1);
Scholl_std=std(tempScholl,1)/sqrt(size(tempScholl,1));
size(Scholl)
resampleX=output.resampleX;

figure;
plot(resampleX,Scholl,'red')
output.meanScholl=Scholl;
output.stdScholl=Scholl_std;
output.tempScholl=tempScholl;
end
