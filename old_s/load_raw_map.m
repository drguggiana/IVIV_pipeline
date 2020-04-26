function load_raw_map(pathName,str,ccr)
%% Load raw map and plot

flipo=str(ccr).sliceOri;

%pathName='C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\L23\Paper\Cellsfor_iviv_paper\Cellsfor_iviv_paper\170613\SW0002\map02'
list=dir([char(pathName) filesep '*.xsg']);
j=1;
load([char(pathName) filesep list(j).name],'-mat');
sr = header.ephys.ephys.sampleRate;%check sample rate
srF = 1/(1000/sr);
samples_per_sweep = header.ephys.ephys.traceLength*sr;
timebase=1/sr:1/sr:samples_per_sweep/sr; %TR2019: timebase
traces=data.ephys.trace_1;%raw ephys trace
%Filter
cutoff      = 1000;     % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
order       = 4;        % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
type        = 'Butter';
traces = lowpassfilt(traces, order, cutoff, sr, type);
ind_traces=reshape(traces,[length(traces)/256 256]);
%Baseline subtraction
base_start=1;
base_end=100;
baseline=ind_traces(base_start*srF:base_end*srF,:);
bs_traces=ind_traces-mean(baseline);%subtract baseline
flipFlag = 0;
flipx = 1;
array=bs_traces;
traceLength = size(array,1 );
mapPattern=header.mapper.mapper.mapPatternArray;
if flipx ==1
    mapPattern=fliplr(mapPattern)
end
 for n=1:numel(mapPattern)
 newArray(:,find(mapPattern==n)) = array(:,n);
 end
startTime = 1;
stopTime = traceLength;
showStart = 99*srF;
showStop = 700*srF;
array = newArray(startTime:stopTime, :);
[rows,cols] = size(array);
totalTime = (rows-1)/sr; 
xTimeAxis = linspace(0, totalTime, rows)';
traceAxis = ( 1 : cols );
[sizeX, sizeY] = size(mapPattern);
yFactor = 500; % offset, in pA
%xFactor = 500
scaleBarTxt = 'pA';
if flipFlag == 1
yFactor = -yFactor; 
end
offsetVector = yFactor * ( 0 : cols-1 );
offsetArray = repmat(offsetVector, rows, 1);
array = array-offsetArray;      
% set up the figure -------------------------------------------------------------
x = .14; 
y = .14; 
w = .8; 
h = .8; 
F=figure;
subplotRows = 1; subplotCols = sizeY; plotnum = 0;
for N = 1:sizeY
startInd = (N-1)*sizeX + 1;
endInd = N*sizeX;
plotnum = plotnum+1;

    %     hsub(plotnum) = subplot(subplotRows,subplotCols,plotnum);
    
%     pos1 = 0.025 + (N - 1)*(0.96/sizeY);
%      pos2 = 0.02;
%     pos3 = 0.05;
%     pos4 = 0.96;

   pos1 = 0.025 + (N - 1)*(0.96/sizeY);
    pos2 = 0.02;
    pos3 = 0.038;
    pos4 = 0.96

    hsub(N) = axes('Position', [pos1 pos2 pos3 pos4]);

%     set(gca, 'Position', );

    plot(xTimeAxis(showStart:showStop), array(showStart:showStop,startInd:endInd),'Color','k');
    hold on;
%                 y1=get(gca,'ylim');
%                 x1= redpeak_start/(srF*100);
%                 hold on;
%                 p1=plot([x1 x1],y1,'--','Color','r');
%                 p1.Color(4) = 0.3;
%                 hold on;
% %                 y1=get(gca,'ylim');
% %                 x1=redpeak_end/(srF*100);
% %                 hold on;
% %                 p2=plot([x1 x1],y1,'--','Color','r');
% %                 p2.Color(4) = 0.3;
% %                 hold on;
%                 %blue vertical lines
%                 y1=get(gca,'ylim');
%                 x1=bluepeak_start/(srF*100);
%                 hold on;
%                 p3=plot([x1 x1],y1,'--','Color','b');
%                 p3.Color(4) = 0.3;
%                 hold on;
%                 y1=get(gca,'ylim');
%                 x1=bluepeak_end/(srF*100);
%                 hold on;
%                 p4=plot([x1 x1],y1,'--','Color','b');
%                 p4.Color(4) = 0.3;
    
minval = min(mean(array(1:100,startInd:endInd)));
maxval = max(mean(array(1:100,startInd:endInd)));
tweakFactor = abs(maxval - minval)*0.05;
yrange = [minval-tweakFactor maxval+tweakFactor];
set(gca, 'YLim', yrange);
set(gca, 'XLim', [(showStart-200)/sr (showStop+200)/sr]);
xlabel('Seconds');
set(gcf,'color','w');
%  if flipo==0
%      set(gca, 'XDir','reverse')
% end   
                
end
          
set(findobj(gcf, 'Type', 'axes'), 'Visible', 'off');
    
% scalebar lines
    Y = mean(array(:,end))+yFactor/4;
% Y = min(get(gca, 'YLim'));
    hscalebar = line([.1 .2], [Y Y]);
    set(hscalebar, 'Color', 'k', 'Tag', 'scaleBarLines');
    hscalebar = line([.1 .1], [Y Y+yFactor/2]);
    set(hscalebar, 'Color', 'k', 'Tag', 'scaleBarLines');

% scalebar text
    ht(1) = text(.12, Y+yFactor/6, '100 ms'); 
    ht(2) = text(.12, Y+yFactor/3, [num2str(yFactor/2) ' ' scaleBarTxt]); 
    set(ht, 'Color', 'k', 'FontSize', 8, 'Tag', 'scaleBarText');
%axis square;
end