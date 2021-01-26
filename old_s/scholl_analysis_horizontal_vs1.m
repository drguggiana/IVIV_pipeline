function output=scholl_analysis_horizontal_vs1(data,sampling,extend)
%VS 07/01/21
%input:
%data - matrix with start and endpoint of each dendritic segment in each row
%sampling - sampling of x coordinate
%extend - interval to be covered
%output.resampleScholl - distribution of crossings at sampling points of x
%output.resampleX - vector of sampling points of x

%remove segments with same x coord
flag=data(:,1)-data(:,2);
tempIndex=find(flag~=0);
data=data(tempIndex,:);

%make absolute of colum one entries smaller than of column 2
tempComp=abs(data(:,1))>abs(data(:,2));
tempIndex=find(tempComp==1);
tempData=data;
data(tempIndex,1)=tempData(tempIndex,2);
data(tempIndex,2)=tempData(tempIndex,1);

flag=data(:,:)./abs(data(:,:));
flag(isnan(flag(:)))=0;
flag=flag(:,1)+flag(:,2);

%right side (x>0)
tempIndex=find(flag>0);
rightCoords=data(tempIndex,:);

%left side (x<0)
tempIndex=find(flag<0);
leftCoords=data(tempIndex,:);

%´segments crossing midline
tempIndex=find(flag==0);
crossingCoords=data(tempIndex,:);
%left side of crossing segments
tempIndex=find(crossingCoords(:)<0);
tempCoords(:,2)=crossingCoords(tempIndex);
tempCoords(:,1)=0;
leftCoords=[leftCoords; tempCoords];
%right side of crossing segments
tempIndex=find(crossingCoords(:)>0);
tempCoords(:,2)=crossingCoords(tempIndex);
tempCoords(:,1)=0;
rightCoords=[rightCoords; tempCoords];

%left scholl distribution
leftX=-leftCoords(:);
leftNums=[ones(size(leftCoords,1),1); -ones(size(leftCoords,1),1)];
[leftX,I]=sort(leftX,'ascend');
leftNums=leftNums(I);
leftScholl=cumsum(leftNums);
%right scholl distribution
rightX=rightCoords(:);
rightNums=[ones(size(rightCoords,1),1); -ones(size(rightCoords,1),1)];
[rightX,I]=sort(rightX,'ascend');
rightNums=rightNums(I);
rightScholl=cumsum(rightNums);
%complete Sholl
completeScholltemp=[leftScholl(end:-1:1);rightScholl];
completeXtemp=[-leftX(end:-1:1);rightX];

%remove multiple entries
uniqueX=unique(completeXtemp);
for i=1:size(uniqueX,1)
    tempIndex=find(completeXtemp==uniqueX(i));
    if uniqueX(i)<0
        completeScholl(i)=completeScholltemp(tempIndex(1));
    else
        completeScholl(i)=completeScholltemp(tempIndex(end));
    end
    completeX(i)=uniqueX(i);
end

%resample
resampleX=[-extend:sampling:extend];
resampleXtemp=[-extend:sampling:-sampling];
for i=1:size(resampleXtemp,2)
    tempIndex=find(completeX>=resampleXtemp(i));
    resampleScholl_left(i)=completeScholl(min(tempIndex));
end
resampleXtemp=[sampling:sampling:extend];
for i=1:size(resampleXtemp,2)
    tempIndex=find(completeX<=resampleXtemp(i));
    resampleScholl_right(i)=completeScholl(max(tempIndex));
end
tempIndex=find(completeX==0);
resampleScholl=[resampleScholl_left completeScholl(tempIndex) resampleScholl_right];

% % for illustration:
% figure;hold
% % plot(rightX,rightScholl)
% plot(completeX,completeScholl)
% plot(resampleX,resampleScholl,'red')

output.resampleScholl=resampleScholl;
output.resampleX=resampleX;
end
