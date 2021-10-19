%% plot visit in cues over time

black_visitOverTime= visitOverTime;
figure()
index= 1;
plot(black_visitOverTime(index).nearTestAIR)
lg=[{black_visitOverTime(index).date}];
hold on
for index= 2:size(black_visitOverTime,2)
    plot(black_visitOverTime(index).nearTestAIR)
    lg= horzcat(lg, {black_visitOverTime(index).date});
end
hold off
legend(lg);

black_visitOverTime= visitOverTime;
totalExp= size(black_visitOverTime,2);
figure();
for index=1:totalExp
    subplot(3,4,index);
    plot(black_visitOverTime(index).nearTestAIR)
    hold on
    plot(black_visitOverTime(index).nearBaseAIR)
    hold off
    lg= [{'Near Test Cues'}, {'Near Base Cue'}];
    legend(lg);
    title(strcat(black_visitOverTime(index).date,' Visits near cues - AIR'));
end

figure();
for index=1:totalExp
    subplot(3,4,index);
    plot(black_visitOverTime(index).nearTestCO2)
    hold on
    plot(black_visitOverTime(index).nearBaseCO2)
    hold off
    lg= [{'Near Test Cues'}, {'Near Base Cue'}];
    legend(lg);
    title(strcat(black_visitOverTime(index).date,' Visits near cues - CO2'));
end

figure();
for index=1:totalExp
    subplot(3,4,index);
    plot(black_visitOverTime(index).nearTestPostCO2)
    hold on
    plot(black_visitOverTime(index).nearBasePostCO2)
    hold off
    lg= [{'Near Test Cues'}, {'Near Base Cue'}];
    legend(lg);
    title(strcat(black_visitOverTime(index).date,' Visits near cues - PostCO2'));
end


%test to plot the data smoothed
plot(smooth([visitOverTime(index).nearTestAIR, visitOverTime(index).nearTestCO2, visitOverTime(index).nearTestPostCO2],0.25), 'red')
hold on
plot(smooth([visitOverTime(index).nearTestAIR, visitOverTime(index).nearTestCO2, visitOverTime(index).nearTestPostCO2],0.95), 'black')
%plot(smooth([visitOverTime(index).nearTestAIR, visitOverTime(index).nearTestCO2, visitOverTime(index).nearTestPostCO2],0.1,'rloess'), 'blue')
plot(smooth([visitOverTime(index).nearTestAIR, visitOverTime(index).nearTestCO2, visitOverTime(index).nearTestPostCO2],0.1), 'blue')
hold off
lg=[{'smooth 1'}, {'smooth 2'}, {'smooth 3'}];
legend(lg);




% plot the counts in Test + base cues together and in plot is AIR + CO2 + postCO2
color= 'black';
clear dataset
dataset=black_visitOverTime;
for index= 1:size(dataset,2)
    dataToPlot= dataset(index).nearTestAIR + dataset(index).nearBaseAIR;
    dataToPlot= horzcat(dataToPlot,  dataset(index).nearTestCO2 + dataset(index).nearBaseCO2);
    dataToPlot= horzcat(dataToPlot,  dataset(index).nearTestPostCO2 + dataset(index).nearBasePostCO2);
    maxValue= max(dataToPlot);
    %plot(dataToPlot, 'color','black')
%     plot(smooth(dataToPlot,0.1), 'color','black')
%     figure()
%     plot(smooth(dataToPlot,0.25), 'color','black')
%     figure()
%     plot(smooth(dataToPlot,0.05, 'loess'), 'color','black')
%     figure()
%     plot(smooth(dataToPlot,0.1, 'loess'), 'color','black')
%     figure()
    plot(smooth(dataToPlot,0.1, 'loess'), 'color','black')
    hold on
    plot([61 61], [0 maxValue/2], ':r')
    plot([121 121], [0, maxValue/2], ':r')
    hold off
    t= char(dataset(index).date);
    title(strcat(t(1:8),'-',t(10:15)));
    outputPath= 'D:\DataFlydra\MosquitoProject\OUTPUT_Graphs\';
    saveas(gcf,strcat(outputPath, 'smooth_',color,'_',mat2str(index),'_countsOverTime.pdf'));
end


color= 'gray40';
clear dataset
dataset=gray40_visitOverTime;
for index= 1:size(dataset,2)
    dataToPlot= dataset(index).nearTestAIR + dataset(index).nearBaseAIR;
    dataToPlot= horzcat(dataToPlot,  dataset(index).nearTestCO2 + dataset(index).nearBaseCO2);
    dataToPlot= horzcat(dataToPlot,  dataset(index).nearTestPostCO2 + dataset(index).nearBasePostCO2);
    maxValue= max(dataToPlot);
    plot(dataToPlot, 'color',[0.6 0.6 0.6])
    hold on
    plot([61 61], [0 maxValue/2], ':r')
    plot([121 121], [0, maxValue/2], ':r')
    hold off
    t= char(dataset(index).date);
    title(strcat(t(1:8),'-',t(10:15)));
    outputPath= 'D:\DataFlydra\MosquitoProject\OUTPUT_Graphs\';
    saveas(gcf,strcat(outputPath, color,'_',mat2str(index),'_countsOverTime.pdf'));
end


color= 'gwt1';
clear dataset
dataset=gwt1_visitOverTime;
for index= 1:size(dataset,2)
    dataToPlot= dataset(index).nearTestAIR + dataset(index).nearBaseAIR;
    dataToPlot= horzcat(dataToPlot,  dataset(index).nearTestCO2 + dataset(index).nearBaseCO2);
    dataToPlot= horzcat(dataToPlot,  dataset(index).nearTestPostCO2 + dataset(index).nearBasePostCO2);
    maxValue= max(dataToPlot);
    plot(dataToPlot, 'color',[0 1 0]*0.5)
    hold on
    plot([61 61], [0 maxValue/2], ':r')
    plot([121 121], [0, maxValue/2], ':r')
    hold off
    t= char(dataset(index).date);
    title(strcat(t(1:8),'-',t(10:15)));
    outputPath= 'D:\DataFlydra\MosquitoProject\OUTPUT_Graphs\';
    saveas(gcf,strcat(outputPath, color,'_',mat2str(index),'_countsOverTime.pdf'));
end


color= 'gct1';
clear dataset
dataset=gct1_visitOverTime;
for index= 1:size(dataset,2)
    dataToPlot= dataset(index).nearTestAIR + dataset(index).nearBaseAIR;
    dataToPlot= horzcat(dataToPlot,  dataset(index).nearTestCO2 + dataset(index).nearBaseCO2);
    dataToPlot= horzcat(dataToPlot,  dataset(index).nearTestPostCO2 + dataset(index).nearBasePostCO2);
    maxValue= max(dataToPlot);
    plot(dataToPlot, 'color',[0 1 0]*0.8)
    hold on
    plot([61 61], [0 maxValue/2], ':r')
    plot([121 121], [0, maxValue/2], ':r')
    hold off
    t= char(dataset(index).date);
    title(strcat(t(1:8),'-',t(10:15)));
    outputPath= 'D:\DataFlydra\MosquitoProject\OUTPUT_Graphs\';
    saveas(gcf,strcat(outputPath, color,'_',mat2str(index),'_countsOverTime.pdf'));
end


color= 'rhue';
clear dataset
dataset=rhue_visitOverTime;
for index= 1:size(dataset,2)
    dataToPlot= dataset(index).nearTestAIR + dataset(index).nearBaseAIR;
    dataToPlot= horzcat(dataToPlot,  dataset(index).nearTestCO2 + dataset(index).nearBaseCO2);
    dataToPlot= horzcat(dataToPlot,  dataset(index).nearTestPostCO2 + dataset(index).nearBasePostCO2);
    maxValue= max(dataToPlot);
    plot(dataToPlot, 'color','red')
    hold on
    plot([61 61], [0 maxValue/2], ':r')
    plot([121 121], [0, maxValue/2], ':r')
    hold off
    t= char(dataset(index).date);
    title(strcat(t(1:8),'-',t(10:15)));
    outputPath= 'D:\DataFlydra\MosquitoProject\OUTPUT_Graphs\';
    saveas(gcf,strcat(outputPath, color,'_',mat2str(index),'_countsOverTime.pdf'));
end


%% Calculate % of trajectories visited ANY CUE
% This section uses data stored in the smryCounts table
% (load_sumary_report) and expData.listTotalIDsNearCue(AIR/CO2)
tmpVisitAIR=zeros(length(expData),1);
tmpVisitCO2=zeros(length(expData),1);
for i= 1:length(expData)
    expName= expData(i).name;
    tmpIndex= find(strcmp(smryTrajectories.date, expData(i).name));
    %calcualte % of trajectories near cue
    tmpVisitAIR(i)= length(expData(i).listTotalIDsNearCueAIR)*100/smryTrajectories.numOfTrajAIR(tmpIndex);
    tmpVisitCO2(i)= length(expData(i).listTotalIDsNearCueCO2)*100/smryTrajectories.numOfTrajCO2(tmpIndex);
end


%% TEST 
c1 = animatedline('Color',[0 .7 .7]);
c2 = animatedline('Color',[0 .05 .05]);



sz= 10000;
x = linspace(0,20,sz);
axis([0 20 -1 1])
for k = 1:length(x)
    % first line
    xk = x(k);
    ysin = sin(xk);
    addpoints(c1,xk,ysin);
    
    % second line
    ycos = cos(xk);
    addpoints(c2,xk,ycos);
    
    % update screen
    drawnow %limitrate
end


cList=[];
ids= [471; 345 ];
for i= ids'
    eval(['c' num2str(i) '=animatedline("LineStyle", "none","Marker", "." )']);
    cList= vertcat(cList, sprintf('c%d',i));
end
sz= 10000;
rndTs= sort(randi([1,sz],[sz,1]));
x = linspace(0,20,sz);
d1= [ids(1)*ones(sz,1), rndTs, sin(x)']; 
rndTs= sort(randi([1,sz],[sz,1]));
d2= [ids(2)*ones(sz,1), rndTs, cos(x)'];
cMatrix= vertcat(d1,d2);
cMatrix= sortrows(cMatrix, 2);

axis([0 20 -1 1])

for k = 1:length(x)
    % first line
    xk = x(k);
    %ysin = sin(xk);
    cIndex= find(str2num(cList(:,2:end)) == cMatrix(k,1));
    cSelected= evalin('base', sprintf(cList(cIndex,:)));
    addpoints(cSelected,x(k),cMatrix(k,3));

    % update screen
    drawnow %limitrate
end



%another way for multiple animated lines (updating them manually)
axis([0 20 -1 1])
for k = 1:length(x)
    % first line
    xk = x(k);
    if (cMatrix(k,1) == 1)
        % first line
        addpoints(c1,xk,cMatrix(k,3));
    else
        % second line
        addpoints(c2,xk,cMatrix(k,3));
    end
    % update screen
    drawnow
end



%% (A.I) LOAD DATA CLEAN
% Load in DataClean only data for any rtajectory longer than flightTimeLimit
% Sort the trajectories points in function of the TS recorded (and move each trajectory at a time
for expIndex= 1:length(dataset)
    dataClean(expIndex).fileName= dataset(expIndex).fileName;
    dataClean(expIndex).type= dataset(expIndex).type;
    dataClean(expIndex).gender= dataset(expIndex).gender;
    dataClean(expIndex).expCues= dataset(expIndex).expCues;
    dataClean(expIndex).attr_id= [];
    dataClean(expIndex).attr_time= [];
    dataClean(expIndex).attr_frame= [];
    dataClean(expIndex).attr_x= [];
    dataClean(expIndex).attr_y= [];
    dataClean(expIndex).attr_z= [];
    dataClean(expIndex).stim= [];
    
    %load only data from trajectories during AIR and CO2
    tempIndexes= find(strcmp(dataset(expIndex).stim(:),'postCO2'));
    closeCO2= min(dataset(expIndex).attr_time(tempIndexes));
    workIndexes= find(dataset(expIndex).attr_time < closeCO2);
    uniqueID= unique(dataset(expIndex).attr_id(workIndexes));
    dataSmry=[];
    for index= 1:length(uniqueID)
        objID= uniqueID(index);
        %Load the frames where appears the current objID
        objTime= find(dataset(expIndex).attr_id(workIndexes) == objID);
        % Estimate the duration of each insect ID trajectory
        duration= get_trajectory_duration(dataset(expIndex).attr_time(workIndexes(objTime)));
        %Consider only if flight duration is bigger than flightTimeLimit seconds
        % and the trajectory was on the lower level of the WT (zThreshold= lim_z/2)
        zThreshold= lim_z/2;
        if duration >= flightTimeLimit && all(dataset(expIndex).attr_z(workIndexes(objTime)) <= zThreshold)
            dataClean(expIndex).attr_id= vertcat(dataClean(expIndex).attr_id, dataset(expIndex).attr_id(workIndexes(objTime)));
            dataClean(expIndex).attr_time= vertcat(dataClean(expIndex).attr_time, dataset(expIndex).attr_time(workIndexes(objTime)));
            dataClean(expIndex).attr_frame= vertcat(dataClean(expIndex).attr_frame, dataset(expIndex).attr_frame(workIndexes(objTime)));
            dataClean(expIndex).attr_x= vertcat(dataClean(expIndex).attr_x, dataset(expIndex).attr_x(workIndexes(objTime)));
            dataClean(expIndex).attr_y= vertcat(dataClean(expIndex).attr_y, dataset(expIndex).attr_y(workIndexes(objTime)));
            dataClean(expIndex).attr_z= vertcat(dataClean(expIndex).attr_z, dataset(expIndex).attr_z(workIndexes(objTime)));
            dataClean(expIndex).stim= vertcat(dataClean(expIndex).stim, dataset(expIndex).stim(workIndexes(objTime)));
        end
    end
end
% ======================

%% (A.II) Pick 100 random trajectories in AIR and 100 random trajectories in CO2 and plot them
% Using dataClean generated above to test this approach

% Select the IDs for both samples (AIR/CO2)
expIndex=1;
sampleSz= 25;

tmpIndexes= find(strcmp(dataClean(expIndex).stim, 'AIR'));
listIDsAIR= unique(dataClean(expIndex).attr_id(tmpIndexes));
sampleIdxs= randsample(1:length(listIDsAIR), sampleSz);
sampleIDsAIR= listIDsAIR(sampleIdxs');

tmpIndexes= find(strcmp(dataClean(expIndex).stim, 'CO2'));
listIDsCO2= unique(dataClean(expIndex).attr_id(tmpIndexes));
sampleIdxs= randsample(1:length(listIDsCO2), sampleSz);
sampleIDsCO2= listIDsCO2(sampleIdxs');

% Group the XZ values for teh selected IDs
trajMatrixAIR=[];
listSizesAIR=[];
for i= 1:length(sampleIDsAIR)
    objID= sampleIDsAIR(i);
    workIndexes= find(dataClean(expIndex).attr_id == objID);
    trajAIR(i).attr_id= objID;
    trajAIR(i).attr_time=  dataClean(expIndex).attr_time(workIndexes);
    trajAIR(i).attr_x=  dataClean(expIndex).attr_x(workIndexes);
    trajAIR(i).attr_y=  dataClean(expIndex).attr_y(workIndexes);
    trajAIR(i).attr_z=  dataClean(expIndex).attr_z(workIndexes);
    listSizesAIR= vertcat(listSizesAIR, [objID, length(workIndexes)]);
    % Group everything as a matrix (to be able to sort by attr_time column)
    trajMatrixAIR= vertcat(trajMatrixAIR, [trajAIR(i).attr_id*ones(listSizesAIR(i,2),1), trajAIR(i).attr_time, trajAIR(i).attr_x,trajAIR(i).attr_y, trajAIR(i).attr_z]);
end

trajMatrixCO2=[];
listSizesCO2=[];
for i= 1:length(sampleIDsCO2)
    objID= sampleIDsCO2(i);
    workIndexes= find(dataClean(expIndex).attr_id == objID);
    trajCO2(i).attr_id= objID;
    trajCO2(i).attr_time=  dataClean(expIndex).attr_time(workIndexes);
    trajCO2(i).attr_x=  dataClean(expIndex).attr_x(workIndexes);
    trajCO2(i).attr_y=  dataClean(expIndex).attr_y(workIndexes);
    trajCO2(i).attr_z=  dataClean(expIndex).attr_z(workIndexes);
    listSizesCO2= vertcat(listSizesCO2, [objID, length(workIndexes)]);
    % Group everything as a matrix (to be able to sort by attr_time column)
    trajMatrixCO2= vertcat(trajMatrixCO2, [trajCO2(i).attr_id*ones(listSizesCO2(i,2),1), trajCO2(i).attr_time, trajCO2(i).attr_x,trajCO2(i).attr_y, trajCO2(i).attr_z]);
end


%GSCATTER FOR AIR
figure()
subplot(2,1,1);
axis([-lim_x lim_x -lim_y lim_y])
gscatter(trajMatrixAIR(:,3), trajMatrixAIR(:,4), trajMatrixAIR(:,1));
hold on
plot(-0.43, -0.1, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 10); 
plot(-0.43, 0.1, 'o', 'MarkerFaceColor',[0.8, 0.8, 0.8], 'MarkerSize', 6);
plot([-0.455, -0.455],[-0.3, 0.3], 'b-');
plot([-0.405, -0.405],[-0.3, 0.3], 'b-');
title('GSCATTER for AIR')
hold off
subplot(2,1,2);
axis([-lim_x lim_x 0 lim_z])
gscatter(trajMatrixAIR(:,3), trajMatrixAIR(:,5), trajMatrixAIR(:,1));
hold on
plot(-0.43, 0, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 10); 
plot([-0.455, -0.455],[0, 0.3], 'b-');
plot([-0.405, -0.405],[0, 0.3], 'b-');
hold off

%GSCATTER FOR CO2
figure()
subplot(2,1,1);
axis([-lim_x lim_x -lim_y lim_y])
gscatter(trajMatrixCO2(:,3), trajMatrixCO2(:,4), trajMatrixCO2(:,1));
hold on
plot(-0.43, -0.1, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 10); 
plot(-0.43, 0.1, 'o', 'MarkerFaceColor',[0.8, 0.8, 0.8], 'MarkerSize', 6);
plot([-0.455, -0.455],[-0.3, 0.3], 'b-');
plot([-0.405, -0.405],[-0.3, 0.3], 'b-');
title('GSCATTER for CO2')
hold off
subplot(2,1,2);
axis([-lim_x lim_x 0 lim_z])
gscatter(trajMatrixCO2(:,3), trajMatrixCO2(:,5), trajMatrixCO2(:,1));
hold on
plot(-0.43, 0, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 10); 
plot([-0.455, -0.455],[0, 0.3], 'b-');
plot([-0.405, -0.405],[0, 0.3], 'b-');
hold off


% Sort both matrix by the timestamp value
%sortedTrajMatrixAIR= sortrows(trajMatrixAIR, 2);
% 
% % Initialize video
% myVideo = VideoWriter(strcat(outputPath,outputFolder,'trajectoriesWithAIR_60fps_Sorted.mp4')); %open video file
% myVideo.FrameRate = 60;  %can adjust this, 5 - 10 works well for me
% open(myVideo)
% 
% % II. plot all the points as animated lines
% axis([-lim_x lim_x 0 lim_z])
% for k = 1:length(sortedTrajMatrixAIR(:,1))
%     % Select the correct animatedLine (id_XXXX) for the given traj ID XXXX
%     %cIndex= find(str2num(idsPlotList(:,4:end)) == trajMatrixAIR(k,1));
%     %currentID= evalin('base', sprintf(idsPlotList(cIndex,:)));   
%     cIndex= find(strcmp(idsPlotList, sprintf('id_%d',sortedTrajMatrixAIR(k,1))));
%     currentID= evalin('base', char(idsPlotList(cIndex)));
%     
%     addpoints(currentID,sortedTrajMatrixAIR(k,3),sortedTrajMatrixAIR(k,4));
%     % update screen
%     drawnow
%     
%     %pause(0.01) %Pause and grab frame
%     frame = getframe(gcf); %get frame
%     writeVideo(myVideo, frame);  
% end
% 
% close(myVideo)
% 
%  clear id_* idsPlotList myVideo

%  %+++++++++++++++++++++++++
% % Different try
% % Find longest trajectory
% maxTrajSz= max(listSizesAIR);
% axis([-lim_x lim_x 0 lim_z])
% % plot trajectories all starting at the same time
% for k = 1:maxTrajSz
%     for i=1:sampleSz
%         %look fo each traj ID
%         id= sampleIDsAIR(i);
%         currentID= evalin('base', sprintf('id_%d',id)); 
%         if k <= listSizesAIR(i)
%             %load its current animatedLine and add points
%             addpoints(currentID,trajAIR(i).attr_x(k),trajAIR(i).attr_z(k));
%         else
%             clearpoints(currentID);
%          end
%     end
%     % once we have load the currents points fo all traj IDS
%     drawnow
% end
clear id_* trajAIR_2 
%  %+++++++++++++++++++++++++
 %+++++++++++++++++++++++++
% animate trjectories AIR
% I. Generate the variables for each ID as animatedLine

 % Find longest trajectory
maxTrajSz= max(listSizesAIR(:,2));
%extraEmptyStps= fix((maxTrajSz - listSizesAIR)/2);
extraEmptyStps= (maxTrajSz - listSizesAIR(:,2));
 
for i= 1: length(listSizesAIR(:,1))   
    eval(['id_' num2str(listSizesAIR(i,1)) '=animatedline("Color", [0.5 0.5 0.5],"LineStyle", "none","Marker", ".")']);
    % add the extra empty space to its attr_x and attr_z attributes
    %select a random starting time (in the animation)
    if listSizesAIR(i,2) == maxTrajSz
        startingPnt=1;
    else
        startingPnt= randperm(extraEmptyStps(i),1);
    end
    trajAIR_2(i).id= sampleIDsAIR(i);
    trajAIR_2(i).attr_time= [zeros(startingPnt,1); trajAIR(i).attr_time; zeros((maxTrajSz-(startingPnt+listSizesAIR(i,2))),1)];
    trajAIR_2(i).attr_x= [zeros(startingPnt,1); trajAIR(i).attr_x; zeros((maxTrajSz-(startingPnt+listSizesAIR(i,2))),1)];
    trajAIR_2(i).attr_y= [zeros(startingPnt,1); trajAIR(i).attr_y; zeros((maxTrajSz-(startingPnt+listSizesAIR(i,2))),1)];
    trajAIR_2(i).attr_z= [zeros(startingPnt,1); trajAIR(i).attr_z; zeros((maxTrajSz-(startingPnt+listSizesAIR(i,2))),1)];
end 


% Initialize video
myVideo = VideoWriter(strcat(outputPath,outputFolder,'trajectoriesWithAIR_60fps_rndSt_2.mp4')); %open video file
myVideo.FrameRate = 60;  %can adjust this, 5 - 10 works well for me
open(myVideo)


axis([-lim_x lim_x 0 lim_z])
% plot trajectories all starting at the same time
hold on
% Plot Odor Source
plot(-0.76, 0.20, 'ro', 'MarkerFaceColor',[0.4940 0.1840 0.5560], 'MarkerSize', 10);
% Plot Visual Cue
plot(-0.43, 0, 'o', 'MarkerFaceColor','r', 'MarkerSize', 10); 
% Plot Volume
plot([-0.500, -0.500],[0, 0.04], 'b-');
plot([-0.360, -0.360],[0, 0.04], 'b-');
plot([-0.500, -0.360],[0.04, 0.04], 'b-');
for k = 1:maxTrajSz
    for i=1:sampleSz
        %look for each traj ID
        currentID= evalin('base', sprintf('id_%d',listSizesAIR(i,1))); 
        if k < listSizesAIR(i,2)
            %load its current animatedLine and add points
            addpoints(currentID,trajAIR_2(i).attr_x(k),trajAIR_2(i).attr_z(k));
        else
            % Change color to indicate that the traj has ended
            currentID.Color= [0.8, 0.8, 0.8];
            %clearpoints(currentID);
         end
    end
    
    
    % once we have load the currents points fo all traj IDS
    drawnow
    
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);  
end
hold off
close(myVideo);


clear id_*

maxTrajSz= max(listSizesCO2(:,2));
extraEmptyStps= (maxTrajSz - listSizesCO2(:,2));
 
for i= 1: length(listSizesCO2(:,1))   
    eval(['id_' num2str(listSizesCO2(i,1)) '=animatedline("Color", [0.5 0.5 0.5], "LineStyle", "none","Marker", ".")']);
    % add the extra empty space to its attr_x and attr_z attributes
    %select a random starting time (in the animation)
    if listSizesCO2(i,2) == maxTrajSz
        startingPnt=1;
    else
        startingPnt= randperm(extraEmptyStps(i),1);
    end
    trajCO2_2(i).id= sampleIDsCO2(i);
    trajCO2_2(i).attr_time= [zeros(startingPnt,1); trajCO2(i).attr_time; zeros((maxTrajSz-(startingPnt+listSizesCO2(i,2))),1)];
    trajCO2_2(i).attr_x= [zeros(startingPnt,1); trajCO2(i).attr_x; zeros((maxTrajSz-(startingPnt+listSizesCO2(i,2))),1)];
    trajCO2_2(i).attr_y= [zeros(startingPnt,1); trajCO2(i).attr_y; zeros((maxTrajSz-(startingPnt+listSizesCO2(i,2))),1)];
    trajCO2_2(i).attr_z= [zeros(startingPnt,1); trajCO2(i).attr_z; zeros((maxTrajSz-(startingPnt+listSizesCO2(i,2))),1)];
end 


% Initialize video
myVideo = VideoWriter(strcat(outputPath,outputFolder,'trajectoriesWithCO2_60fps_rndSt_3.mp4')); %open video file
myVideo.FrameRate = 60;  %can adjust this, 5 - 10 works well for me
open(myVideo)


axis([-lim_x lim_x 0 lim_z])
% plot trajectories all starting at the same time
hold on
% Plot Odor Source
plot(-0.76, 0.20, 'ro', 'MarkerFaceColor',[0.4940 0.1840 0.5560], 'MarkerSize', 10);
% Plot Visual Cue
plot(-0.43, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 10); 
% Plot Volume
plot([-0.500, -0.500],[0, 0.04], 'b-');
plot([-0.360, -0.360],[0, 0.04], 'b-');
plot([-0.500, -0.360],[0.04, 0.04], 'b-');
for k = 1:maxTrajSz
    for i=1:sampleSz
        %look fo each traj ID
        currentID= evalin('base', sprintf('id_%d',listSizesCO2(i,1))); 
        if k < listSizesCO2(i,2)
            %load its current animatedLine and add points
            addpoints(currentID,trajCO2_2(i).attr_x(k),trajCO2_2(i).attr_z(k));
        else
            % Change color to indicate that the traj has ended
            currentID.Color= [0.8, 0.8, 0.8];         end
    end
    % once we have load the currents points fo all traj IDS
    drawnow
    
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);  
end
hold off
    
close(myVideo);






%+++++++++++++++++++++++++
 
% animate trjectories CO2
% I. Generate the variables for each ID as animatedLine
idsPlotList=[];
for id= sampleIDsCO2'
    eval(['id_' num2str(id) '=animatedline("LineStyle", "-","Marker", ".")']);
    idsPlotList= vertcat(idsPlotList, {sprintf('id_%d',id)});
end

% defaults for the marker
hAx=gca;      % get the axis handle
sz=10;        % size of marker
clr='b';      % color
hS=scatter(hAx,nan,nan,sz,clr);  % initial point won't show but creates handle


% Initialize video
myVideo = VideoWriter(strcat(outputPath,outputFolder,'trajectoriesWithCO2_60fps_2.mp4')); %open video file
myVideo.FrameRate = 60;  %can adjust this, 5 - 10 works well for me
open(myVideo)

% II. plot all the points as animated lines
axis([-lim_x lim_x 0 lim_z])
for k = 1:length(trajMatrixCO2(:,1))
    % Select the correct animatedLine (id_XXXX) for the given traj ID XXXX
    currentID= evalin('base', sprintf('id_%d',trajMatrixCO2(k,1)));
    %cIndex= find(str2num(idsPlotList(:,4:end)) == trajMatrixAIR(k,1));
    %currentID= evalin('base', sprintf(idsPlotList(cIndex,:)));
    addpoints(currentID,trajMatrixCO2(k,3),trajMatrixCO2(k,4));
    %Update the marker Position
    set(hS,'xdata',trajMatrixCO2(k,3),'ydata',trajMatrixCO2(k,4))      
    % update screen
    drawnow
    
    %pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);  
    
end

close(myVideo)

 clear id_*

 
% =======================================
% =======================================

%% load trajectories behavior when the CO2 is released (from 5 sec before CO2 release to 5 sec after CO2 release)
expIndex=1;
timeDuration= 1000; %1 second
tempIndexes= find(strcmp(dataClean(expIndex).stim(:),'CO2'));
releaseCO2= min(dataClean(expIndex).attr_time(tempIndexes));
timeIndexes= find(dataClean(expIndex).attr_time >= (releaseCO2 -  timeDuration*5) & dataClean(expIndex).attr_time <= (releaseCO2 + timeDuration*5)); 

dataInTimePeriod.attr_id= dataClean(expIndex).attr_id(timeIndexes);
dataInTimePeriod.attr_time= dataClean(expIndex).attr_time(timeIndexes);
dataInTimePeriod.attr_frame= dataClean(expIndex).attr_frame(timeIndexes);
dataInTimePeriod.attr_x= dataClean(expIndex).attr_x(timeIndexes);
dataInTimePeriod.attr_y= dataClean(expIndex).attr_y(timeIndexes);
dataInTimePeriod.attr_z= dataClean(expIndex).attr_z(timeIndexes);
dataInTimePeriod.attr_stim= dataClean(expIndex).stim(timeIndexes);

% Select X trajectories randomly
sampleSz= 100;
listIDs= unique(dataInTimePeriod.attr_id);
sampleIdx= randsample(1:length(listIDs), sampleSz);
sampleIDs= listIDs(sampleIdx);


trajMatrix=[];
listSizes=[];
for i= 1:length(sampleIDs)
    objID= sampleIDs(i);
    workIndexes= find(dataInTimePeriod.attr_id == objID);
    traj(i).attr_id= objID;
    traj(i).attr_time=  dataInTimePeriod.attr_time(workIndexes);
    traj(i).attr_x=  dataInTimePeriod.attr_x(workIndexes);
    traj(i).attr_z=  dataInTimePeriod.attr_z(workIndexes);
    listSizes= vertcat(listSizes, length(workIndexes));
    % Group everything as a matrix (to be able to sort by attr_time column)
    trajMatrix= vertcat(trajMatrix, [traj(i).attr_id*ones(listSizes(i),1), traj(i).attr_time, traj(i).attr_x, traj(i).attr_z]);
end

% Sort by attr_time
trajTimeSorted= sortrows(trajMatrix, 2);

% Plot all trajectories together
gscatter(trajMatrix(:,3), trajMatrix(:,4), trajMatrix(:,1))
figure()
gscatter(trajTimeSorted(:,3), trajTimeSorted(:,4), trajTimeSorted(:,1))

% = PLot over time =
% I. Generate the variables for each ID as animatedLine
idsPlotList=[];
for id= sampleIDs'
    eval(['id_' num2str(id) '=animatedline("LineStyle", "none","Marker", ".")']);
    idsPlotList= vertcat(idsPlotList, {sprintf('id_%d',id)});
end

% II. plot all the points as animated lines
axis([-lim_x lim_x 0 lim_z])
for k = 1:length(trajTimeSorted(:,1))
    % Select the correct animatedLine (id_XXXX) for the given traj ID XXXX 
    %cIndex= find(str2num(idsPlotList(:,4:end)) == trajTimeSorted(k,1));
    %currentID= evalin('base', sprintf(idsPlotList(cIndex,:)));
    cIndex= find(strcmp(idsPlotList, sprintf('id_%d',trajTimeSorted(k,1))));
    currentID= evalin('base', char(idsPlotList(cIndex)));
    addpoints(currentID,trajTimeSorted(k,3),trajTimeSorted(k,4));
    % update screen
    drawnow
end


%% PLot traj in a given odor, when the trajectories have left the CO2 plume area
% Load in DataClean only data for any trajectory longer than flightTimeLimit
% Check the # of points for each trajectory and match the # of points (with
% value X, Z value= 999) with the biggest trajectory size
for expIndex=1:length(dataset) 
    expIndex=1;
    dataClean(expIndex).fileName= dataset(expIndex).fileName;
    dataClean(expIndex).type= dataset(expIndex).type;
    dataClean(expIndex).gender= dataset(expIndex).gender;
    dataClean(expIndex).expCues= dataset(expIndex).expCues;
    dataClean(expIndex).attr_id= [];      
    dataClean(expIndex).attr_time= [];
    dataClean(expIndex).attr_frame= [];
    dataClean(expIndex).attr_x= [];
    dataClean(expIndex).attr_y= [];
    dataClean(expIndex).attr_z= [];
    dataClean(expIndex).stim= [];
    
    %load only data from trajectories during AIR and CO2
    tempIndexes= find(strcmp(dataset(expIndex).stim(:),'postCO2'));
    closeCO2= min(dataset(expIndex).attr_time(tempIndexes)); 
    workIndexes= find(dataset(expIndex).attr_time < closeCO2);
    uniqueID= unique(dataset(expIndex).attr_id(workIndexes));
    dataSmry=[]
    for index= 1:length(uniqueID)
        objID= uniqueID(index);
        %Load the frames where appears the current objID
        objTime= find(dataset(expIndex).attr_id(workIndexes) == objID);
        % Estimate the duration of each insect ID trajectory
        duration= get_trajectory_duration(dataset(expIndex).attr_time(workIndexes(objTime)));
        %Consider only if flight duration is bigger than flightTimeLimit seconds
        % and the trajectory was on the lower level of the WT (zThreshold= lim_z/2)
        zThreshold= lim_z/2;
        if duration >= flightTimeLimit && all(dataset(expIndex).attr_z(workIndexes(objTime)) <= zThreshold)
            dataClean(expIndex).attr_id= vertcat(dataClean(expIndex).attr_id, dataset(expIndex).attr_id(workIndexes(objTime)));      
            dataClean(expIndex).attr_time= vertcat(dataClean(expIndex).attr_time, dataset(expIndex).attr_time(workIndexes(objTime)));
            dataClean(expIndex).attr_frame= vertcat(dataClean(expIndex).attr_frame, dataset(expIndex).attr_frame(workIndexes(objTime)));
            dataClean(expIndex).attr_x= vertcat(dataClean(expIndex).attr_x, dataset(expIndex).attr_x(workIndexes(objTime)));
            dataClean(expIndex).attr_y= vertcat(dataClean(expIndex).attr_y, dataset(expIndex).attr_y(workIndexes(objTime)));
            dataClean(expIndex).attr_z= vertcat(dataClean(expIndex).attr_z, dataset(expIndex).attr_z(workIndexes(objTime)));
            dataClean(expIndex).stim= vertcat(dataClean(expIndex).stim, dataset(expIndex).stim(workIndexes(objTime)));           
        end
    end
end   

    expIndex= 1;
    timeDuration= 1000; %1 second
    tempIndexes= find(strcmp(dataClean(expIndex).stim(:),'CO2'));
    releaseCO2= min(dataClean(expIndex).attr_time(tempIndexes));
    timeIndexes= find(dataClean(expIndex).attr_time >= (releaseCO2 - timeDuration) & dataClean(expIndex).attr_time <= (releaseCO2 + timeDuration*2)); 
    %[a,b]= groupcounts(dataset(expIndex).stim(timeIndexes))    %test
    
    dataInTimePeriod.attr_id= dataClean(expIndex).attr_id(timeIndexes);
    dataInTimePeriod.attr_time= dataClean(expIndex).attr_time(timeIndexes);
    dataInTimePeriod.attr_frame= dataClean(expIndex).attr_frame(timeIndexes);
    dataInTimePeriod.attr_x= dataClean(expIndex).attr_x(timeIndexes);
    dataInTimePeriod.attr_y= dataClean(expIndex).attr_y(timeIndexes);
    dataInTimePeriod.attr_z= dataClean(expIndex).attr_z(timeIndexes);
    dataInTimePeriod.attr_stim= dataClean(expIndex).stim(timeIndexes);
    
    % gplotmatrix(dataInTimePeriod.attr_x, dataInTimePeriod.attr_z, dataInTimePeriod.attr_id)
    
    % Pick a set of Trajectories ID
    listIDs= [3652, 3654, 3822, 4293, 4305, 3616, 3647, 3696, 4439]; 
   
    listSizes= [];
    for i= 1:length(listIDs)
        objID= listIDs(i);
        workIndexes= find(dataInTimePeriod.attr_id == objID);
        traj(i).attr_id= objID;
        traj(i).attr_x=  dataInTimePeriod.attr_x(workIndexes);
        traj(i).attr_z=  dataInTimePeriod.attr_z(workIndexes);
        listSizes= vertcat(listSizes, length(workIndexes));
    end
    maxSz= max(listSizes);
    for i= 1:length(listIDs)
       if listSizes(i) < maxSz
           extraSz= maxSz - listSizes(i);
           traj(i).attr_x= vertcat(traj(i).attr_x, 999*ones(extraSz,1));
           traj(i).attr_z= vertcat(traj(i).attr_z, 999*ones(extraSz,1));
       end
    end
    a1 = animatedline('Color',[.1 .1 .1], 'Marker','o', 'LineStyle', 'none');
    a2 = animatedline('Color',[.4 .4 .4], 'Marker','o', 'LineStyle', 'none');    
    a3 = animatedline('Color',[.6 .6 .6], 'Marker','o', 'LineStyle', 'none');
    a4 = animatedline('Color',[.8 .1 .8], 'Marker','o', 'LineStyle', 'none');
    durationInFrames= linspace(1,maxSz, maxSz);
    %durationInFrames= linspace(1,fps*3, fps*3);

    xlim([-lim_x, lim_x]);
    ylim([0, lim_z]);
    for k= 1:length(durationInFrames)   %index= 1:length(listIDs)
        addpoints(a1,traj(3).attr_x(k),traj(3).attr_z(k));
        addpoints(a2,traj(7).attr_x(k),traj(7).attr_z(k));
        addpoints(a3,traj(9).attr_x(k),traj(9).attr_z(k));
        addpoints(a4,traj(5).attr_x(k),traj(5).attr_z(k));
        pause(0.12)
        %update Screen
        drawnow limitrate
    end



% ====================================

%% PLot trajs 1 sec before and 2 second after the CO2 release
%for expIndex= 1
ts_startAIR= 1584472390; 
ts_startCO2= 1584475991;
ts_endCO2= 1584483191;
ts_endAIR= 1584490391; 
%---------------

for expIndex=1:length(dataset)
    expIndex= 1;
    
    dataClean(expIndex).fileName= dataset(expIndex).fileName;
    dataClean(expIndex).type= dataset(expIndex).type;
    dataClean(expIndex).gender= dataset(expIndex).gender;
    dataClean(expIndex).expCues= dataset(expIndex).expCues;
    dataClean(expIndex).attr_id= [];      
    dataClean(expIndex).attr_time= [];
    dataClean(expIndex).attr_frame= [];
    dataClean(expIndex).attr_x= [];
    dataClean(expIndex).attr_y= [];
    dataClean(expIndex).attr_z= [];
    dataClean(expIndex).stim= [];
    
    %load only data from trajectories during AIR and CO2
    tempIndexes= find(strcmp(dataset(expIndex).stim(:),'postCO2'));
    closeCO2= min(dataset(expIndex).attr_time(tempIndexes)); 
    workIndexes= find(dataset(expIndex).attr_time < closeCO2);
    uniqueID= unique(dataset(expIndex).attr_id(workIndexes));
    dataSmry=[]
    for index= 1:length(uniqueID)
        objID= uniqueID(index);
        %Load the frames where appears the current objID
        objTime= find(dataset(expIndex).attr_id(workIndexes) == objID);
        % Estimate the duration of each insect ID trajectory
        duration= get_trajectory_duration(dataset(expIndex).attr_time(workIndexes(objTime)));
        %Consider only if flight duration is bigger than flightTimeLimit seconds
        if duration >= flightTimeLimit
            dataClean(expIndex).attr_id= vertcat(dataClean(expIndex).attr_id, dataset(expIndex).attr_id(workIndexes(objTime)));      
            dataClean(expIndex).attr_time= vertcat(dataClean(expIndex).attr_time, dataset(expIndex).attr_time(workIndexes(objTime)));
            dataClean(expIndex).attr_frame= vertcat(dataClean(expIndex).attr_frame, dataset(expIndex).attr_frame(workIndexes(objTime)));
            dataClean(expIndex).attr_x= vertcat(dataClean(expIndex).attr_x, dataset(expIndex).attr_x(workIndexes(objTime)));
            dataClean(expIndex).attr_y= vertcat(dataClean(expIndex).attr_y, dataset(expIndex).attr_y(workIndexes(objTime)));
            dataClean(expIndex).attr_z= vertcat(dataClean(expIndex).attr_z, dataset(expIndex).attr_z(workIndexes(objTime)));
            dataClean(expIndex).stim= vertcat(dataClean(expIndex).stim, dataset(expIndex).stim(workIndexes(objTime)));
            
        end
    end
end   
    %[~,idx]= sort(dataSmry(:,2));
    %dataSmrySorted= dataSmry(idx,:);
    %Obtain moment the CO2 was released, and take the indexes around that timestamp
    expIndex= 1;
    timeDuration= 1000; %1 second
    tempIndexes= find(strcmp(dataClean(expIndex).stim(:),'CO2'));
    releaseCO2= min(dataClean(expIndex).attr_time(tempIndexes));
    timeIndexes= find(dataClean(expIndex).attr_time >= (releaseCO2 - timeDuration) & dataClean(expIndex).attr_time <= (releaseCO2 + timeDuration*2)); 
    %[a,b]= groupcounts(dataset(expIndex).stim(timeIndexes))    %test
    
    dataInTimePeriod.attr_id= dataClean(expIndex).attr_id(timeIndexes);
    dataInTimePeriod.attr_time= dataClean(expIndex).attr_time(timeIndexes);
    dataInTimePeriod.attr_frame= dataClean(expIndex).attr_frame(timeIndexes);
    dataInTimePeriod.attr_x= dataClean(expIndex).attr_x(timeIndexes);
    dataInTimePeriod.attr_y= dataClean(expIndex).attr_y(timeIndexes);
    dataInTimePeriod.attr_z= dataClean(expIndex).attr_z(timeIndexes);
    dataInTimePeriod.attr_stim= dataClean(expIndex).stim(timeIndexes);
    
    uniqueIDs= unique(dataInTimePeriod.attr_id);    
    for index= 1:length(uniqueIDs)
        h= figure();
        objID= uniqueIDs(index);
        %Load the frames where appears the current objID
        workIndexes= find(dataInTimePeriod.attr_id == objID);
        % Sort the indexes for the objID in funtion of the timestamp
        %[out,idx] = sort(dataClean(expIndex).attr_time(workIndexes));
        if  any(dataInTimePeriod.attr_time(workIndexes) <= releaseCO2) && ~any(dataInTimePeriod.attr_time(workIndexes) >= releaseCO2)            
            plot(dataInTimePeriod.attr_x(workIndexes), dataInTimePeriod.attr_z(workIndexes), 'b');
        elseif ~any(dataInTimePeriod.attr_time(workIndexes) <= releaseCO2) && any(dataInTimePeriod.attr_time(workIndexes) >= releaseCO2)            
            plot(dataInTimePeriod.attr_x(workIndexes), dataInTimePeriod.attr_z(workIndexes), 'g');
        end
        hold on
        plot(dataInTimePeriod.attr_x(workIndexes(1)), dataInTimePeriod.attr_z(workIndexes(1)), '+b');
        plot(dataInTimePeriod.attr_x(workIndexes(end)),dataInTimePeriod.attr_z(workIndexes(end)), 'or');

        if any(dataInTimePeriod.attr_time(workIndexes) <= releaseCO2) && any(dataInTimePeriod.attr_time(workIndexes) >= releaseCO2)
            disp(strcat('ID- ',mat2str(objID),'-was detected in AIR and CO2'));
            changeOdor= find(dataInTimePeriod.attr_time(workIndexes) <= releaseCO2);
            changeOdor= changeOdor(end);
            plot(dataInTimePeriod.attr_x(workIndexes(changeOdor)), dataInTimePeriod.attr_z(workIndexes(changeOdor)), '*k');
        end
        
        hold off
        title(strcat(mat2str(objID),'Trajectory'));
        xlim([-lim_x, lim_x]);
        ylim([0, lim_z]);
        uiwait(h);
    end
    

% 
%     uniqueIDs= unique(dataInTimePeriod.attr_id);    
%     for index= 1:length(uniqueIDs)
%         h= figure();
%         objID= uniqueIDs(index);
%         %Load the frames where appears the current objID
%         workIndexes= find(dataClean(expIndex).attr_id(timeIndexes) == objID);
%         % Sort the indexes for the objID in funtion of the timestamp
%         %[out,idx] = sort(dataClean(expIndex).attr_time(workIndexes));
%         if  any(dataClean(expIndex).attr_time(timeIndexes(workIndexes)) <= releaseCO2) && ~any(dataClean(expIndex).attr_time(timeIndexes(workIndexes)) >= releaseCO2)            
%             plot(dataClean(expIndex).attr_x(workIndexes), dataClean(expIndex).attr_z(workIndexes), 'b');
%         elseif ~any(dataClean(expIndex).attr_time(timeIndexes(workIndexes)) <= releaseCO2) && any(dataClean(expIndex).attr_time(workIndexes) >= releaseCO2)            
%             plot(dataClean(expIndex).attr_x(workIndexes), dataClean(expIndex).attr_z(workIndexes), 'g');
%         end
%         hold on
%         plot(dataClean(expIndex).attr_x(workIndexes(1)), dataClean(expIndex).attr_z(workIndexes(1)), '+b');
%         plot(dataClean(expIndex).attr_x(workIndexes(end)), dataClean(expIndex).attr_z(workIndexes(end)), 'or');
% 
%         if any(dataClean(expIndex).attr_time(workIndexes) <= releaseCO2) && any(dataClean(expIndex).attr_time(workIndexes) >= releaseCO2)
%             disp(strcat('ID- ',mat2str(objID),'-was detected in AIR and CO2'));
%             changeOdor= find(dataClean(expIndex).attr_time(workIndexes) <= releaseCO2);
%             changeOdor= changeOdor(end);
%             plot(dataClean(expIndex).attr_x(workIndexes(changeOdor)), dataClean(expIndex).attr_z(workIndexes(changeOdor)), '*k');
%         end
%         
%         hold off
%         uiwait(h);
%     end
%     



%% Identify ID that are in AIR AND CO2 parts (All traj plotted at the same time without animation)
for expIndex=1:length(dataset)
    expIndex= 10;
    indexWithCO2 = find(strcmp(dataset(expIndex).stim(:),'CO2'));
    releaseCO2= min(dataset(expIndex).attr_time(indexWithCO2));

    indexPrevCO2= find(dataset(expIndex).attr_time >= (releaseCO2 - timeDuration) & dataset(expIndex).attr_time < releaseCO2);
    indexWithCO2= find(dataset(expIndex).attr_time >= releaseCO2 & dataset(expIndex).attr_time <= (releaseCO2 + timeDuration)); 
    timeIndexes= vertcat(indexPrevCO2, indexWithCO2);
    uniqueIDsAIR= unique(dataset(expIndex).attr_id(indexPrevCO2)); 
    uniqueIDsCO2= unique(dataset(expIndex).attr_id(indexWithCO2));
    % compare the unique ID values in both groups
    IDsInBoth= ismember(uniqueIDsCO2,uniqueIDsAIR);
    disp(expIndex)
    %groupcounts(IDsInBoth)
    IDsToPlot= uniqueIDsCO2(IDsInBoth);
    
    %now obtain their X-Z values
    id= IDsToPlot(1);

    h=figure()
    hold on
    for i= 1:length(IDsToPlot)
        traj(i).id= IDsToPlot(i);
        trajIndexes= find(dataset(expIndex).attr_id(timeIndexes) == traj(i).id);
        traj(i).tList= dataset(expIndex).attr_time(timeIndexes(trajIndexes));
        traj(i).xList= dataset(expIndex).attr_x(timeIndexes(trajIndexes));
        traj(i).zList= dataset(expIndex).attr_z(timeIndexes(trajIndexes));
        plot(traj(i).xList,traj(i).zList);
        plot(traj(i).xList(1),traj(i).zList(1),'+');
        plot(traj(i).xList(end),traj(i).zList(end), 'o');
        
    end
    xlim([-0.1, 0.5]);
    ylim([0, 0.3]);
    hold off 
    for index=1:length(xList)
        h=plot(xList(index), zList(index), '*r');
        if tList(index)>= releaseCO2
            plot(0.15,0.15, '+b')
        end
        pause(0.12)
        delete(h)
    end
    hold off
    
    clear indexPrevCO2 indexWithCO2 IDsInBoth
    
end

% ====================================
%=====================================

%% Plot a given trajectory ID
% Black vs white - 20200317_121244
% [14941,14978,15332,16374,16464,17072,17246,17795,18186]
%Green vs White – '20200620_121920'
% 9150	9295	9331	9478	9488	9537	9547
trajID= 9331; 
indexesObjID= find(dataset(15).attr_id == trajID);
objXYZ= [dataset(15).attr_x(indexesObjID), dataset(15).attr_y(indexesObjID), dataset(15).attr_z(indexesObjID)];
expCues= dataset(15).expCues;
expCues(2,1)={'green'};
plot_trajectory_2D_v3(trajID, objXYZ, 'black', expCues, 1.6, 0)


%% Generate your own speed calculator

for odorChecked=[{'AIR'},{'CO2'},{'postCO2'}]
    for expIndex=1:length(dataset)
        odorIndex= find(strcmp(dataset(expIndex).stim, odorChecked));
        data=[dataset(expIndex).attr_id(odorIndex), dataset(expIndex).attr_time(odorIndex), dataset(expIndex).attr_x(odorIndex), dataset(expIndex).attr_y(odorIndex), dataset(expIndex).attr_z(odorIndex)];
        %select all the IDs from the vector
        uniqueID=unique(data(:,1));
        %Transpose from a column matrix to a row matrix
        uniqueID=uniqueID';
        % Initialice local variables
        totalTraj=0; 
        flightTimeList= nan(length(uniqueID),1);
        totalFlightTime=0;
        speedList=[];

        % For each insect ID, check if 
        for index= 1:length(uniqueID)
            objID= uniqueID(index);
            %Load the frames where appears the current objID
            objTime= find(data(:,1) == objID);
            % Estimate the duration of each insect ID trajectory
            duration= get_trajectory_duration(data(objTime(:),2));

            %Consider only if flight duration is bigger than flightTimeLimit seconds
            if duration >= flightTimeLimit
                % Calculate the lenght traveled by trajectory ID in each axis
                distX= sum(abs(diff(data(objTime(:),3))));
                distY= sum(abs(diff(data(objTime(:),4))));
                distZ= sum(abs(diff(data(objTime(:),5))));

                % Estimate the mean velocity
                distT= sqrt(distX^2 + distY^2 + distZ^2);
                speed= distT/duration;

                speedList= [speedList, speed]; 
                
            end
        end
        
        disp(strcat('-- ', int2str(expIndex), ' ---- ', odorChecked))
        disp(strcat('max speed:', sprintf('%.6f',max(speedList))))
        if strcmp(odorChecked, 'AIR')
            meanSpeedList(expIndex, 1)= mean(speedList);
        elseif strcmp(odorChecked, 'CO2')
            meanSpeedList(expIndex, 2)= mean(speedList);
        elseif strcmp(odorChecked, 'postCO2')
            meanSpeedList(expIndex, 3)= mean(speedList);
        end
    end
end 

plot(meanSpeedList)
lg=[{'AIR'}, {'CO2'}, {'postCO2'}]
legend(lg, 'Location', 'northeast');




%% Load Speed values XYZ and estimate average flightspeed (per trajectory)
% I dont know which type of speed is generated by FLYDRA and stored in the
% .h5 file as xvel, yvel, zvel. 

inputPath= strcat(inputPath,subFolder);
cd(inputPath)
filesList=dir('*mainbrain.h5');

% Move to the MATLAB workspace
cd(workspace);

%initialize list for tested color used and position ob tyhe base color cue in WT
baseColorIndexList=zeros(length(filesList),1);
testedColorList= cell(length(filesList),1);

loadFullDataset=true; % Variable used to test with partial data
fileNameList={};
% ===========================  
% Load and clean the raw data from Flydra
for expIndex= 1:length(filesList)
        
        % Load file name to work with 
        %fileName= '20200223_110547.mainbrain.h5';
        fileName= filesList(expIndex).name;
        filePath= strcat(inputPath,fileName);
        fileNameList=vertcat(fileNameList, {fileName(1:15)});
        %disp (filePath);
        % Load all the information from the FLYDRA .h5 file    
        %[attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z]= load_data_from_file(filePath, loadFullDataset);
        [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z, attr_xvel, attr_yvel, attr_zvel]= load_PosAndSpeed_from_file(filePath, loadFullDataset);
        %tempData= [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z];
        
        % If using .JSON files
        % data=jsondecode(fileread('data.txt'));
        
        if contains(subFolder, 'AnStephensi')
            % Load experiment settings for the an.Stephensi mosquitoes
            [cuesSetup,ts_startAIR,ts_startCO2,ts_endCO2,ts_endAIR, mType, mGender]=  load_exp_settings_anStephensi(fileName);
        elseif contains(subFolder, 'CxQuinquefasciatus')
            [cuesSetup,ts_startAIR,ts_startCO2,ts_endCO2,ts_endAIR, mType, mGender]=  load_exp_settings_cxQuinquefasciatus(fileName);
        else
            % Load the experiment settings (cues positions & exp timestamps)
            [cuesSetup,ts_startAIR,ts_startCO2,ts_endCO2,ts_endAIR, mType, mGender]=  load_exp_settings(fileName);            
        end
        % --- Clean the data loaded from the H5 file ---
        % Erase the X,Y and Z values that are out of the test section limits
        %[attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z]= erase_pos_outside_WT_v2(attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z, lim_x, lim_y, lim_z);
        [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z, attr_xvel, attr_yvel, attr_zvel]= erase_pos_outside_WT_v3(attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z, attr_xvel, attr_yvel, attr_zvel, lim_x, lim_y, lim_z);
        
        % Erase the entries recorded before and after the mass flow
        % controller script is launched and the entries with timeStamp == 0.0.
        expIndexes= find(attr_time(:) >= ts_startAIR & attr_time(:) <= ts_endAIR & attr_time(:) ~=0.0);
        % Erase all insectIDs that have been detected in only 1 frame (or
        % below 0.x seconds (is it worth it to keep it?)
        % --- ---

        %Fill the dataset with the data from current experiment
        dataset(expIndex).fileName=    fileName;
        dataset(expIndex).type=        mType;
        dataset(expIndex).gender=       mGender;
        dataset(expIndex).expCues=    cuesSetup;
        
        dataset(expIndex).attr_id=     attr_id(expIndexes);
        dataset(expIndex).attr_time=   attr_time(expIndexes);
        dataset(expIndex).attr_frame=  attr_frame(expIndexes);
        dataset(expIndex).attr_x=      attr_x(expIndexes);
        dataset(expIndex).attr_y=      attr_y(expIndexes);
        dataset(expIndex).attr_z=      attr_z(expIndexes);
        dataset(expIndex).attr_xvel=      attr_xvel(expIndexes);
        dataset(expIndex).attr_yvel=      attr_yvel(expIndexes);
        dataset(expIndex).attr_zvel=      attr_zvel(expIndexes);

        % Estimate the timestamps values for odor stimulus ON and OFF
        % attr_time contains timestamp epochs (in seconds)
        indexPrevCO2= find(dataset(expIndex).attr_time(:) < ts_startCO2);
        indexWithCO2= find(dataset(expIndex).attr_time(:) >= ts_startCO2 & dataset(expIndex).attr_time(:) < ts_endCO2-1);
        indexPostCO2= find(dataset(expIndex).attr_time(:) >= ts_endCO2);
                
        dataset(expIndex).stim(indexPrevCO2,1)= {'AIR'};
        dataset(expIndex).stim(indexWithCO2,1)= {'CO2'};
        dataset(expIndex).stim(indexPostCO2,1)= {'postCO2'};
              
        if length(cuesSetup(:,1)) == 3
            %Load the position where the baseCue is and the color tested in experiment
            baseCueIndex=find(strcmp(cuesSetup([2,3],1), baseColor));
            baseColorIndexList(expIndex)= baseCueIndex;
            testedColorList(expIndex)= cuesSetup(find(~strcmp(cuesSetup([2,3],1), baseColor))+1,1);
        end
end


%Clear temporary variables from workspace
clear fileName loadFullDataset indexPrevCO2 indexWithCO2 indexPostCO2
clear attr_id attr_time attr_frame attr_x attr_y attr_z
clear cuesSetup ts_startAIR ts_startCO2 ts_endCO2 ts_endAIR mType mGender
clear filesList expIndexes


%%
%% load trajectories inside volume over time
dataTest= [ 3	4838	1584476039.01424;
            3	4838	1584476039.03286;
            3	4838	1584476039.04423;
            1	4866	1584476058.05937;
            1	4866	1584476058.07798;
            1	4866	1584476058.09337;
            1	4866	1584476058.11117;
            2	12116	1584480894.97694;
            2	12116   1584480894.99161;
            2	12116	1584480895.06194;
            1	12140	1584480911.8959;
            1	12140	1584480911.91058;
            2	12140	1584480907.52329;
            2	12140	1584480907.53741;
            1	12240	1584480992.85194;
            1	12240	1584480992.87036;
            1	12240	1584480992.88505;
            1	12240	1584480992.90977;
            2	12240	1584480989.12926;
            2	12240	1584480989.15273;
            2	12240	1584480989.16331;
            3	12271	1584481034.97117;
            3	12271	1584481034.97978;
            1	12279	1584481056.71152]
initialTS= 1584475991.00703;
finalTS= 1584483188.72342;
numGrps= 6;

% % Initialize the counters
% p1 = zeros(length(filesList),numGrps); 
% p2 = zeros(length(filesList),numGrps);
% p3 = zeros(length(filesList),numGrps); 
% p4 = zeros(length(filesList),numGrps);

for fileIndex= 1:length(filesList)
    fileName= filesList(fileIndex).name;
    disp(strcat(' - Working with file: ', {' '}, fileName));
    filesName(fileIndex)= {fileName(1:8)};
    % data loaded structure: [vCuePosition	objID	timeStamp NaN initialTimeStamp finalTimeStamp]
    dataFromExcel = table2array(readtable(strcat(filesPath, fileName)));
    if isempty(dataFromExcel)
        disp(strcat('   -> No counts detected in file: ',{' '}, fileName));
        % Add zeros to the counts inside each time group 
        p1(fileIndex,:)= 0;
        p2(fileIndex,:)= 0;
        p3(fileIndex,:)= 0;
        p4(fileIndex,:)= 0;
    else
            %Load initial and final timestamps for the odor used
            initialTS= dataFromExcel(1,4);    
            finalTS= dataFromExcel(1,5);
            %disp(strcat('  - data from excel:',num2str(dataFromExcel(10,3))));
            %disp(strcat('  - data from excel 2:',num2str(dataFromExcel(1430,3))));
            grpSz= (finalTS - initialTS)/numGrps;
            grpThreshold= zeros(1, numGrps+1);
            for i=0:numGrps
                grpThreshold(i+1)= initialTS+(grpSz*i);
                %grpThreshold(i+1)= (grpSz*i);

            end
            %g=groupcounts(dataFromExcel(:,3), initialTS+grpThreshold );
            for grpIndex= 1:numGrps %length(grpThreshold)
                %disp(strcat('  - grpIndex:',num2str(grpIndex)));
                %disp(strcat('  - threshold:', num2str(grpThreshold(grpIndex))));

                %Select the indexes from the data that verify the time group
                %conditions
                dataIndexes= find(dataFromExcel(:,3) < grpThreshold(grpIndex+1) & dataFromExcel(:,3) >= grpThreshold(grpIndex));
                p1(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 1);
                p2(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 2);
                p3(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 3);
                p4(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 4);
            end
    end
end






%% estimate amount of trajectories in a random volume in the wind tunnel
% volume size
radius= 0.07;
h= 0.04;

% Generate a random volume position per experiment
for expIndex=1: length(dataset)
    %Pick a random spot inside the XY plane of the test section (the volume
    %must be fully inside the test section
    cuesX= -lim_x + cell2mat(dataset(expIndex).expCues(2,2));
    tmpLimX= lim_x - radius;
    tmpLimY= lim_y - radius;
    tmpX = (cuesX+radius) + rand*(tmpLimX-(cuesX+radius));
    tmpY = -tmpLimY + rand*(tmpLimY-(-tmpLimY));

    centerList(expIndex,:)= [tmpX tmpY];   
end


odorChecked= 'AIR';
for expIndex= 1:length(dataset)
    fileName= strcat(dataset(expIndex).fileName(1:15),'_randomVol_',odorChecked,'.xls');
    outputFile= strcat(outputPath,outputFolder,'analysisData_randomVolume\', fileName);
    % Load exp Data related to odorChecked value
    indexWithOdor = find(strcmp(dataset(expIndex).stim(:), odorChecked));
    dataEntry.expCues= dataset(expIndex).expCues;
    dataEntry.attr_id= dataset(expIndex).attr_id(indexWithOdor);
    dataEntry.attr_time= dataset(expIndex).attr_time(indexWithOdor);
    dataEntry.attr_x= dataset(expIndex).attr_x(indexWithOdor);
    dataEntry.attr_y= dataset(expIndex).attr_y(indexWithOdor);
    dataEntry.attr_z= dataset(expIndex).attr_z(indexWithOdor);

    % Table containing all the counts information. 
    % 1st row is filled with 0s and will be erased at the end of this function 
    % [positionX  OBJID, instantTimeStamp]
    insectCtrInPos= zeros(1,3);
    allCheckXYZ= 0;

    % Find intial and final timestamps associated to this dataEntry
    initialTS= min(dataEntry.attr_time);
    finalTS= max(dataEntry.attr_time);

    %Pick the trajectory of a given insect
    for objID= unique(dataEntry.attr_id)'
        %Load the frames where appears the current objID
        objFrames= find(dataEntry.attr_id == objID);
        % Align the new objFrames indexes to their real position in dataEntry
        %Load the frames where appears the current objID
        %calculate the duration of the flight
        duration= get_trajectory_duration(dataEntry.attr_time(objFrames));

        if duration >= flightTimeLimit
            % For each visual cue used, find its [X,Y] position and the number
            % of times the insect objID has been near the cue
            center=centerList(expIndex,:);
            % check if the insect is inside the given volume for each of
            % the axis separately
            checkX= find((dataEntry.attr_x(objFrames) > (center(1) - radius)) & (dataEntry.attr_x(objFrames) < (center(1) + radius)));
            checkY= find((dataEntry.attr_y(objFrames) > (center(2) - radius)) & (dataEntry.attr_y(objFrames) < (center(2) + radius)));
            checkZ= find(dataEntry.attr_z(objFrames) < h);
            % Then compare the indexes to see which indexes are inside the
            % volume in the 3 axis at the same time. these indexes are
            % related to the objFrames subset and not the full
            % DataEntry set (*)
            checkXY= intersect(checkX, checkY);
            checkXYZ= intersect(checkXY, checkZ);

            if ~isempty(checkXYZ) 
                % (*) Align CHECK XYZ indexes to real dataEntry indexes
                %indexesInDataEntry= checkXYZ(:) + indexStimStart - 1;
                indexesInDataEntry= objFrames(checkXYZ);
                %Create a row for each of the counts inside the volume
                % [positionX  OBJID, instantTimeStamp]
                %tableRows= [zeros(length(indexesInDataEntry),1), dataEntry.attr_id(indexesInDataEntry), dataEntry.attr_time(indexesInDataEntry)];
                tableRows= [zeros(length(indexesInDataEntry),1), dataEntry.attr_id(indexesInDataEntry), dataEntry.attr_time(indexesInDataEntry)];
                % Add the cue order position to all the entries in
                % tableRows (1-4 for 4 visual cues) 
                tableRows(:,1)= 3;
                %Add the new rows generated to the table
                insectCtrInPos= vertcat(insectCtrInPos, tableRows);
                allCheckXYZ= vertcat(allCheckXYZ, indexWithOdor(indexesInDataEntry));
            end
        end
    end
    %disp(strcat(' - MIN value: ',num2str(min(insectCtrInPos(:)))));
    %disp(strcat(' - MAX value: ',num2str(max(insectCtrInPos(:)))));
    % Erase the first row (it is not real data, just a row full of zeros)
    if sum(insectCtrInPos(1,:)) == 0
        insectCtrInPos(1,:)= [];      
    end
    if allCheckXYZ(1,1) == 0
        allCheckXYZ(1)=[];
    end
    
    % to test
    ctsInVol(expIndex).ids= insectCtrInPos;
    ctsInVol(expIndex).xyz= allCheckXYZ;
    
    
       
    if nnz(insectCtrInPos)
        % Convert the matrix into a table
        insectCtrTable = array2table(insectCtrInPos, 'VariableNames',{'vCuePosition','objID','timeStamp'});
        % Write data in xlsx file
        writetable(insectCtrTable, outputFile, 'Sheet', odorChecked, 'Range', 'A1');
    
        % Write the initial and final timestamps the Excel File
        writetable(table(initialTS, finalTS), outputFile, 'Sheet', odorChecked, 'Range', 'D1')
        % WARNING!: the indexes are related to the subset of data for the odor analyzed (odorChecked). 
        % To test values in full expriment ==> datase do dataset(expIndex).attr_id(indexesOdor(allCheckedXYZ))
        writetable(table(allCheckXYZ), outputFile, 'Sheet', odorChecked, 'Range', 'F1')
        writetable(table(center), outputFile, 'Sheet', odorChecked, 'Range', 'G1')
        % Count the occurences for each visual Cue position by count the number of identical 
        % subscripts (1,2...n) in insectCtrInPos(:,1)
        countsPerPosition= accumarray(insectCtrInPos(:,1),1);
        countsPerPosition= countsPerPosition';

    else
        % If there is no insect inside the volumes
        disp('    -------> Warning! No individuals detected inside the volume for the experiment above ^');
        %Create a file with empty data, will be easier to post analyze if
        %file isempty but created in the folder
        insectCtrTable = array2table(insectCtrInPos, 'VariableNames',{'vCuePosition','objID','timeStamp'});
        % Wrtie data in xlsx file
        writetable(insectCtrTable, outputFile, 'Sheet', odorChecked, 'Range', 'A1');
        writetable(table(center), outputFile, 'Sheet', odorChecked, 'Range', 'G1')

        countsPerPosition= [0, 0];
    end  
end
clear ctsInVol

%%


%% Compare PI in ANY cue vs PI visiting both cues (for 1st hour vs full c02 duration)
% THis is a test using GRAY and data stored in grayExpData_1stHrCO2 and grayExpData_FullCO2
for i= 1:28
    mean1hrCO2_anyCue(i)= mean(expData1stHrCO2(i).PrefIndexTotalIDsCO2);
    mean1hrCO2_bothCues(i)=  mean(expData1stHrCO2(i).PrefIndexBothCuesCO2);
    meanFullCO2_anyCue(i)=  mean(expDataFullCO2(i).PrefIndexTotalIDsCO2);
    meanFullCO2_bothCues(i)=  mean(expDataFullCO2(i).PrefIndexBothCuesCO2); 
end

meanFullCO2_anyCue(:) - mean1hrCO2_anyCue(:)
meanFullCO2_bothCues(:) - mean1hrCO2_bothCues(:)

%data in 1:12 has 2 hours of CO2 and data in 13:28 has 1 hour of CO2
[h p ] = ttest(meanFullCO2_anyCue(1:12), mean1hrCO2_anyCue(1:12))
[h p ] = ttest2(meanFullCO2_anyCue(1:12), mean1hrCO2_anyCue(1:12))

% If p value is close to 0 there is no variation between population
% behavior between taking 2 hours or only the first 1hour.





%%
% Estimating Pref Index per trajectory ID that has visited any cue by
% counting the number of frames where a trajectory ID is near each cue:
% Traj ID has spent some frames only in Base cue:   PI= -1
% Traj ID has spent some frames only in Tested cue: PI= 1
% Traj ID has spent some frames near both cues:     PI= (framesInTest-framesInBase)/(framesInTest+framesInBase)
%warning: It needs the BaseColorIndexList to be generated or loaded from report file 
sheetName= 'exp_with_base_color_white';
[smryTrajectories, smryCounts, smryPIs]= load_summary_report(outputPath, 'Mosquito_Project_Report.xls', sheetName);

% From the summary report file, find the indexes for the experiments we
% want to analyze ('gray2.5', 'gray4.0', 'gray4.5', 'gray6.5', 'gray9.5')
cIndexes= find(strncmp(smryCounts.testedColor, 'gray', 4));
baseColorIndexList= smryCounts.baseColor_Pos(cIndexes);
testColorList= unique(smryCounts.testedColor(cIndexes));

%(1.) Estimate time spent near visual cues for each one of trhe
%trajectories
odorChecked='CO2';
% Load files
%filesPath= strcat(outputPath, outputFolder,'analysisData_oldVersion\');
filesPath= strcat(outputPath, outputFolder,'analysisData_newVersion\');
cd(filesPath);
filesList=dir(strcat('*countsInsideCueVol_',odorChecked,'.xlsx'));
cd(workspace);
for fileIndex= 1:length(filesList)
    fileName= filesList(fileIndex).name;
    disp(strcat(' - Working with file: ', {' '}, fileName));
    filesName(fileIndex)= {fileName(1:8)};
    % data loaded structure: [vCuePosition	objID	timeStamp NaN initialTimeStamp finalTimeStamp]
    dataFromExcel = table2array(readtable(strcat(filesPath, fileName)));
    expData(fileIndex).name= fileName;
    % Initialize the correct fields in function of the odor used
    if strcmp(odorChecked, 'AIR')
        expData(fileIndex).listTotalIDsNearCueAIR=[];
        expData(fileIndex).PrefIndexTotalIDsAIR=[];        
        expData(fileIndex).listIDsInBothCuesAIR=[];
        expData(fileIndex).PrefIndexBothCuesAIR=[];
    elseif strcmp(odorChecked, 'CO2')
        expData(fileIndex).listTotalIDsNearCueCO2=[];
        expData(fileIndex).PrefIndexTotalIDsCO2=[];
        expData(fileIndex).listIDsInBothCuesCO2=[];
        expData(fileIndex).PrefIndexBothCuesCO2=[];
    end
    if isempty(dataFromExcel)
        disp(strcat('   -> No counts detected in file: ',{' '}, fileName));
    else
        %Initialize the object that will contain the information for the current experiment
        tmpPIs=[];
        tmpIDs1=[];
        tmpIDs2=[];
        tmpIDsInBoth=[];
        % To load timestamp data (time spent by insect near the cue
        % Sorted the the data by cue visited ID timestamp
        sortedData= sortrows(dataFromExcel(:,1:3));
        if length(unique(sortedData(:,1)))== 1
            % Find all the Traj IDs in the file
            tmpIDs1= unique(sortedData(:,2));
            if baseColorIndexList(fileIndex) == sortedData(1,1)
                % All the counts where detected near the Base color cue (PI == -1 for all traj IDs)
                tmpPIs= ones(length(tmpIDs),1)*(-1);
            else
                % All the counts where detected near the Tested color cue (PI == 1 for all traj IDs)
                tmpPIs= ones(length(tmpIDs),1);
            end
        else
            % Split the data in function to the position values (1, 2)
            split= find(sortedData(:,1)==2, 1);
            % Count # of frames where each Traj ID visit cue 1 or 2
            [tmpCts1, tmpIDs1]= groupcounts(sortedData(1:(split-1),2));
            [tmpCts2, tmpIDs2]= groupcounts(sortedData(split:end,2));
            
            % Find the traj IDs visiting both cues
            [tmpIDsInBoth, inTmpIDs1, inTmpIDs2]= intersect(tmpIDs1, tmpIDs2);
            
            if any(tmpIDsInBoth)
                % Assign PI values in function of where the Base cue was placed                
                if baseColorIndexList(fileIndex)== 1
                    %Base color was placed in position 1 (+Y axis) in current experiment
                    for i= 1:length(tmpIDsInBoth)
                        tmpPIs(i,1)= (tmpCts2(inTmpIDs2(i))-tmpCts1(inTmpIDs1(i)))/(tmpCts2(inTmpIDs2(i))+tmpCts1(inTmpIDs1(i)));
                    end
                    % Remove ID visiting both cues as they are grouped in tmpIDsInBoth
                    % (Erase duplicate ID when it has visited both cues)
                    tmpIDs1= setdiff(tmpIDs1, tmpIDsInBoth);
                    tmpIDs2= setdiff(tmpIDs2, tmpIDsInBoth);
                    % Assign their Pref Index value
                    tmpPIsOnly1= ones(length(tmpIDs1),1)*(-1);  %PI= -1
                    tmpPIsOnly2= ones(length(tmpIDs2),1);       %PI= 1
                else
                    %Base color was placed in position 2 (-Y axis) in current experiment
                    for i= 1:length(tmpIDsInBoth)
                        tmpPIs(i,1)= (tmpCts1(inTmpIDs1(i))-tmpCts2(inTmpIDs2(i)))/(tmpCts1(inTmpIDs1(i))+tmpCts2(inTmpIDs2(i)));
                    end
                    % Remove ID visiting both cues as they are grouped in tmpIDsInBoth
                    % (Erase duplicate ID when it has visited both cues)
                    tmpIDs2= setdiff(tmpIDs2, tmpIDsInBoth);
                    tmpIDs1= setdiff(tmpIDs1, tmpIDsInBoth);
                    tmpPIsOnly1= ones(length(tmpIDs1),1);       %PI= 1
                    tmpPIsOnly2= ones(length(tmpIDs2),1)*(-1);  %PI= -1                                       
                end
            else
                % ANY trajectory have visited BOTH cues
                if baseColorIndexList(fileIndex)== 1
                    %Base color was placed in position 1 (+Y axis) in current experiment
                    tmpPIsOnly1= ones(length(tmpIDs1),1)*(-1);  %PI= -1
                    tmpPIsOnly2= ones(length(tmpIDs2),1);       %PI= 1
                else
                    %Base color was placed in position 2 (-Y axis) in current experiment
                    tmpPIsOnly1= ones(length(tmpIDs1),1);       %PI= 1
                    tmpPIsOnly2= ones(length(tmpIDs2),1)*(-1);  %PI= -1
                end
            end
        end
        %Add the IDs and the their time spent visiting the cue to the structure
        if strcmp(odorChecked, 'AIR')
            expData(fileIndex).listTotalIDsNearCueAIR=[tmpIDs1; tmpIDs2; tmpIDsInBoth]';
            expData(fileIndex).PrefIndexTotalIDsAIR=[tmpPIsOnly1; tmpPIsOnly2; tmpPIs];
            expData(fileIndex).listIDsInBothCuesAIR=tmpIDsInBoth';
            expData(fileIndex).PrefIndexBothCuesAIR=tmpIDs';
        elseif strcmp(odorChecked, 'CO2')           
            expData(fileIndex).listTotalIDsNearCueCO2=[tmpIDs1; tmpIDs2; tmpIDsInBoth]';
            expData(fileIndex).PrefIndexTotalIDsCO2=[tmpPIsOnly1; tmpPIsOnly2; tmpPIs]';
            expData(fileIndex).listIDsInBothCuesCO2=tmpIDsInBoth';
            expData(fileIndex).PrefIndexBothCuesCO2=tmpPIs';
        end
    end
end
clear tmpIDs1 tmpIDs2 tmpIDsBoth tmpPIsOnly1 tmpPIsOnly2 tmpPIs tmpCts1 tmpCts2
clear tmpIDsInBoth inTmpIDs1 inTmpIDs2 split
 

%Generate a SCATTER PLOT for this type of mutation
figure()
for typeIndex=1:length(testColorList)
    % Find which entries in expData are related to each Tested color
    % expIndex give position regarding the cIndexes used to generate expData
    expIndex= find(strcmp(smryCounts.testedColor(cIndexes), testColorList(typeIndex)));
    
    PIsPerTypeCO2=[expData(expIndex).PrefIndexBothCuesCO2];
    allPIsPerTypeCO2= [expData(expIndex).PrefIndexTotalIDsCO2];
    meanValues(typeIndex)= mean(PIsPerTypeCO2);
    meanAllValues(typeIndex)= mean(allPIsPerTypeCO2);
    subplot(2,3,typeIndex)
    scatter(ones(size(PIsPerTypeCO2,2),1), PIsPerTypeCO2, 10, 'jitter', 'on', 'MarkerFaceColor', [0, 0.5, 0], 'MarkerFaceAlpha', 0.5);
    hold on
    scatter(ones(size(allPIsPerTypeCO2,2),1)*2, allPIsPerTypeCO2, 10, 'jitter', 'on', 'MarkerFaceColor', [0, 0, 0.5], 'MarkerFaceAlpha', 0.5);
    % plot also the means values 
    y= [mean(PIsPerTypeCO2), mean(PIsPerTypeCO2)];
    plot([0.9, 1.1], y, 'g', 'LineWidth',5,'MarkerEdgeColor',[0 0.5 0]);
    y= [mean(allPIsPerTypeCO2), mean(allPIsPerTypeCO2)];
    plot([1.9, 2.1], y, 'b', 'LineWidth',5,'MarkerEdgeColor',[0 0 0.5]);
    hold off
    %lg= [{' With AIR'}, {'With CO2'}];
    lg= [{'BothCues'} {'AnyCue'}];
    legend(lg);
    title(strcat('(Only trajID in both cues) PI near Tested color (',testColorList(typeIndex),')'));
    xlim([0,3]);
    ylim([-1.2, 1.2]);
    ylabel('PI')
end


% ====================================
%=====================================





% ====================================
%=====================================
%========= generate % of time spent near cue per time of trajectory (for each
%trajecotry individually)

%(1.) Estimate time spent near visual cues for each one of trhe
%trajectories
colorToAnalyze= 'white';
baseColorCue= 'white';
TestColorCue= 'black';
odorChecked='CO2';
% Load files
%filesPath= strcat(outputPath, outputFolder,'analysisData_oldVersion\');
filesPath= strcat(outputPath, outputFolder,'analysisData\');
cd(filesPath);
filesList=dir(strcat('*countsInsideCueVol_',odorChecked,'.xlsx'));
cd(workspace);
for fileIndex= 1:length(filesList)       
    fileName= filesList(fileIndex).name;
    disp(strcat(' - Working with file: ', {' '}, fileName));
    filesName(fileIndex)= {fileName(1:8)};
    % data loaded structure: [vCuePosition	objID	timeStamp NaN initialTimeStamp finalTimeStamp]
    dataFromExcel = table2array(readtable(strcat(filesPath, fileName)));
    expData(fileIndex).name= fileName;
    % Initialize the correct fields in function of the odor used
    if strcmp(odorChecked, 'AIR')
        expData(fileIndex).listIDsInCueAIR=[];
        expData(fileIndex).listTimeInCueAIR=[];
        expData(fileIndex).listTrajTimeAIR=[];
        expData(fileIndex).percentNearCueAIR=[];
        expDataAllTrajSz(fileIndex).listIDsInCueAIR=[];
        expDataAllTrajSz(fileIndex).listTimeInCueAIR=[];
        expDataAllTrajSz(fileIndex).listTrajTimeAIR=[];
        expDataAllTrajSz(fileIndex).percentNearCueAIR=[];
    elseif strcmp(odorChecked, 'CO2')
        expData(fileIndex).listIDsInCueCO2=[];
        expData(fileIndex).listTimeInCueCO2=[];
        expData(fileIndex).listTrajTimeCO2=[];
        expData(fileIndex).percentNearCueCO2=[];
        expDataAllTrajSz(fileIndex).listIDsInCueCO2=[];
        expDataAllTrajSz(fileIndex).listTimeInCueCO2=[];
        expDataAllTrajSz(fileIndex).listTrajTimeCO2=[];
        expDataAllTrajSz(fileIndex).percentNearCueCO2=[];
    end
    if isempty(dataFromExcel)
        disp(strcat('   -> No counts detected in file: ',{' '}, fileName));
    else
        %Initialize the object that will contain the information for the current experiment            
        tmpListTime=[];
        tmpIndex= 1;
        % To load timestamp data (time spent by insect near thew cue
        % Sorted the the data by cue visited ID timestamp
        sortedData= sortrows(dataFromExcel(:,1:3));
        if length(unique(sortedData(:,1)))== 1
            %If we have counts in only one position, check if is the position assigned to ColorToAnalyze cue
            % if condition not followed, we analyze all the points in file
            if ~strcmp(colorToAnalyze,baseColorCue)
                % Analyzing for TestedColorCue
                if baseColorIndexList(expIndex) == unique(sortedData(:,1))
                    continue             
                end
            else
                % Analyzing for BaseColorCue
                if baseColorIndexList(expIndex) ~= unique(sortedData(:,1))
                    continue             
                end                
            end
        else%if length(unique(sortedData(:,1)))>= 1
            % Split the data in function to the position values (1, 2)
            split= find(sortedData(:,1)==2, 1);
            %disp(strcat('elseif', num2str(split)));
            % Find the first appearence of the position 2 
            if strcmp(colorToAnalyze,baseColorCue)
                % Analyzing data near Base Color cue
                if baseColorIndexList(fileIndex)== 1
                    %Base color was placed in position 1 (+Y axis) in current experiment
                    sortedData= sortedData(1:(split-1),:);
                else
                    %Base color was placed in position 2 (-Y axis) in current experiment
                    sortedData= sortedData(split:end,:);
                end       
            else
                % Analyzing data near Tested Color Cue
                if baseColorIndexList(fileIndex)== 1
                    %Base color was placed in position 1 (+Y axis) in current experiment
                    sortedData= sortedData(split:end,:);
                else
                    %Base color was placed in position 2 (-Y axis) in current experiment
                    sortedData= sortedData(1:(split-1),:);
                end
            end
        end
        %Estimates the time spend for each ID in the pos where was TESTED cue
        tmpIndex=1;
        for objID= unique(sortedData(:,2))'
            idIndexes= find(sortedData(:,2) ==objID);
            % Calculate the Delta time between timestamps values for the id 
            tsDiff= diff(sortedData(idIndexes,3));
            % find which  sequential counts didn't happen "~=consecutively"
            out= find(tsDiff > 0.05);
            t=0;
            k=1;
            if any(out)
                for i= out
                    % the insect left the volume at a given moment
                    t= t + sum(tsDiff(k:(i-1)));
                    k= i+1;
                end
                t= t + sum(tsDiff(k:end));
            else
                %the insects has been inside the volume all the time
                t= sum(tsDiff);
            end
            tmpListTime(1,tmpIndex)= t;
            tmpIndex= tmpIndex+1;
        end
        %Add the IDs and the their time spent visiting the cue to the structure
        if strcmp(odorChecked, 'AIR')
            expData(fileIndex).listIDsInCueAIR=unique(sortedData(:,2))';
            expData(fileIndex).listTimeInCueAIR=tmpListTime;
            expDataAllTrajSz(fileIndex).listIDsInCueAIR=unique(sortedData(:,2))';
            expDataAllTrajSz(fileIndex).listTimeInCueAIR=tmpListTime;     
        elseif strcmp(odorChecked, 'CO2')
            expData(fileIndex).listIDsInCueCO2=unique(sortedData(:,2))';
            expData(fileIndex).listTimeInCueCO2=tmpListTime;
            expDataAllTrajSz(fileIndex).listIDsInCueCO2=unique(sortedData(:,2))';
            expDataAllTrajSz(fileIndex).listTimeInCueCO2=tmpListTime;
        end
    end
end
clear tempIDs indexes dataFromExcel tmpIndex tmpListTime t k idIndexes tsDiff        
 

%(2.) Once we have generated the time spent near the cues, estimate the
%amount of time spent by each of the trajectories ids
%tic;
odorChecked='CO2';
%expDataAllTrajSz= expData;
for expIndex= 1:length(dataset)
    % Initialize temp var
    tmpIndex=1;
    tmpListTrajTime=[];
    disp(strcat(' - working with: ',expData(expIndex).name));
    %Find the indexes for AIR, CO2 and postCO2 to create a working subdataset
    indexOdor = find(strcmp(dataset(expIndex).stim(:), odorChecked)); 
    tempDataset= [dataset(expIndex).attr_id(indexOdor) dataset(expIndex).attr_time(indexOdor)];
    % Select all the IDs from the dataset
     if strcmp(odorChecked, 'AIR')
         tmpListIDs= expDataAllTrajSz(expIndex).listIDsInCueAIR;
     elseif strcmp(odorChecked, 'CO2')
        tmpListIDs= expDataAllTrajSz(expIndex).listIDsInCueCO2;
     end
    for index= 1:length(tmpListIDs)
        disp(index)
        objID= tmpListIDs(index);
        %Load the frames where appears the current objID
        objTime= find(tempDataset(:,1) == objID);
        if empty(objTime)
            % Remove artifact
            disp(objID)
        end
        % Estimate the duration of each insect ID trajectory
        duration= get_trajectory_duration(tempDataset(objTime,2));
%         %Consider only if flight duration is bigger than flightTimeLimit seconds
%         if (duration >= flightTimeLimit)
%             visitIndex= find(expData(expIndex).listIDsInP1 == objID);
%             if any(visitIndex)
%                 expData(expIndex).totalTimeTrajP1(visitIndex)= duration;
%             end
%             visitIndex= find(expData(expIndex).listIDsInP2 == objID);
%             if any(visitIndex)
%                 expData(expIndex).totalTimeTrajP2(visitIndex)= duration;
%             end
%             expData(expIndex).allTrajID(1,index)=  objID;
%             expData(expIndex).allTrajtotalTime(1,index)=  duration;
%         end
        tmpListTrajTime(1,tmpIndex)= duration;
        tmpIndex= tmpIndex+1;
    end
    if strcmp(odorChecked, 'AIR')
        expDataAllTrajSz(expIndex).listTrajTimeAIR= tmpListTrajTime;
    elseif strcmp(odorChecked, 'CO2')
        expDataAllTrajSz(expIndex).listTrajTimeCO2= tmpListTrajTime; 
    end
end
clear tmpIndex tmpListTrajTime


% Estimate the % of time the traj has spent in the cue per experiment
for expIndex=1:size(expData,2)
    % Group experiments data in function of type of mosquitoes
    timeTrajNearCueAIR=[];
    timeTrajNearCueCO2=[];
    tmpNorm= bsxfun(@rdivide, (expDataAllTrajSz(expIndex).listTimeInCueAIR*100),expDataAllTrajSz(expIndex).listTrajTimeAIR);
    expDataAllTrajSz(expIndex).percentNearCueAIR= tmpNorm;
    tmpNorm= bsxfun(@rdivide, (expDataAllTrajSz(expIndex).listTimeInCueCO2*100),expDataAllTrajSz(expIndex).listTrajTimeCO2);
    expDataAllTrajSz(expIndex).percentNearCueCO2= tmpNorm;
end



%toc;

% Estimate the % of time the traj has spent in the cue
% if strcmp(odorChecked, 'AIR')
%     tmpNorm= bsxfun(@divide, (expData(expIndex).listTimeInCueAIR(:)*100),expData(expIndex).listTrajTimeAIR(:);
% elseif strcmp(odorChecked, 'CO2')
%     tmpNorm= bsxfun(@divide, (expData(expIndex).listTimeInCueCO2(:)*100),expData(expIndex).listTrajTimeCO2(:);
% end 

%mosqType=['m0';'m1';'wt';'l1';'l2';'l3';'l4';'l5';'l6'];
mosqTypes={'m0';'m0';'m1';'wt';'l1';'l4';'l3';'l5';'l2';'l6';'l1';'l6';'l2';'l5';'wt';'l3';'l4';'l5';'l3';'l4';'wt'};
uniqueTypes={'m0';'m1';'wt';'l1';'l2';'l3';'l4';'l5';'l6'};
figure()
for typeIndex=1:size(uniqueTypes,1)
    % Group experiments data in function of type of mosquitoes
    typeValue= uniqueTypes(typeIndex);
    i= find(strcmp(mosqTypes, typeValue));
    timeTrajNearCueAIR=[];
    timeTrajNearCueCO2=[];
    for expIndex=i'
        % FOR AIR
        tmpNorm= bsxfun(@rdivide, (expDataAllTrajSz(expIndex).listTimeInCueAIR*100),expDataAllTrajSz(expIndex).listTrajTimeAIR);
        timeTrajNearCueAIR= [timeTrajNearCueAIR, tmpNorm];
        % Add percentage values to the strucutre
        expDataAllTrajSz(expIndex).percentNearCueAIR= tmpNorm;
        % FOR CO2
        tmpNorm= bsxfun(@rdivide, (expDataAllTrajSz(expIndex).listTimeInCueCO2*100),expDataAllTrajSz(expIndex).listTrajTimeCO2);
        timeTrajNearCueCO2= [timeTrajNearCueCO2, tmpNorm];
        % Add percentage values to the strucutre
        expDataAllTrajSz(expIndex).percentNearCueCO2= tmpNorm;
        
    end

    %Generate a SCATTER PLOT for this type of mutation
    subplot(3,3,typeIndex)
    scatter(ones(size(timeTrajNearCueAIR,2),1), timeTrajNearCueAIR, 10, 'jitter', 'on', 'MarkerFaceColor', [0, 0.5, 0], 'MarkerFaceAlpha', 0.5);
    hold on
    scatter(ones(size(timeTrajNearCueCO2,2),1)*2, timeTrajNearCueCO2, 10, 'jitter', 'on', 'MarkerFaceColor', [0.5, 0, 0], 'MarkerFaceAlpha', 0.5);
    % plot also the means values 
    y= [mean(timeTrajNearCueAIR), mean(timeTrajNearCueAIR)];
    plot([0.5, 1.5], y, 'g', 'LineWidth',5,'MarkerEdgeColor',[0 0.5 0]);
    y= [mean(timeTrajNearCueCO2), mean(timeTrajNearCueCO2)];
    plot([1.5, 2.5], y, 'r', 'LineWidth',5,'MarkerEdgeColor',[0.5 0 0]);
    hold off
    lg= [{' With AIR'}, {'With CO2'}];
    %lg= [{' With AIR'}, {'With CO2'}, {'With AIR (postCO2)'}];
    legend(lg);
    title(strcat('% of trajectory time detected near black cue (',typeValue,')'));
    xlim([0,3]);
    ylim([0, 100]);
    ylabel('% of Trajectories')
    %typeIndex= typeIndex+1;
end    
% ================
% ==================


% DO AS BEFORE BUT ONLY WITH TRAJECTORIES > 1.5SEC
%(2.) Once we have generated the time spent near the cues, estimate the
%amount of time spent by each of the trajectories ids
odorChecked='CO2';
for expIndex= 1:length(dataset)
    % Initialize temp var
    %tmpIndex=1;
    tmpListTrajTime=[];
    
    %Find the indexes for AIR, CO2 and postCO2 to create a working subdataset
    indexOdor = find(strcmp(dataset(expIndex).stim(:), odorChecked)); 
    tempDataset= [dataset(expIndex).attr_id(indexOdor) dataset(expIndex).attr_time(indexOdor)];
    % For each trajectory ID check if it is longer than flightTimeLimit
    % and (if yes) add it to the total for the experiment
     if strcmp(odorChecked, 'AIR')
         tmpListIDs= expData(expIndex).listIDsInCueAIR;
     elseif strcmp(odorChecked, 'CO2')
        tmpListIDs= expData(expIndex).listIDsInCueCO2;
     end
    for tmpIndex= 1:length(tmpListIDs)
        objID= tmpListIDs(tmpIndex);
        %Load the frames where appears the current objID
        objTime= find(tempDataset(:,1) == objID);
        if isempty(objTime)
            % Artifact
            disp(objID);
            duration=0;
        else
            % Estimate the duration of each insect ID trajectory
            duration= get_trajectory_duration(tempDataset(objTime,2));
        end
        %Consider only if flight duration is bigger than flightTimeLimit seconds
        if (duration <= flightTimeLimit)
              duration=NaN;
        end
        tmpListTrajTime(1,tmpIndex)= duration;
        %tmpIndex= tmpIndex+1;
    end
    if strcmp(odorChecked, 'AIR')
        IDsToKeep= ~isnan(tmpListTrajTime);
        expData(expIndex).listIDsInCueAIR= expData(expIndex).listIDsInCueAIR(IDsToKeep);
        expData(expIndex).listTimeInCueAIR= expData(expIndex).listTimeInCueAIR(IDsToKeep);
        expData(expIndex).listTrajTimeAIR= tmpListTrajTime(IDsToKeep);
    elseif strcmp(odorChecked, 'CO2')
        IDsToKeep= ~isnan(tmpListTrajTime);
        expData(expIndex).listIDsInCueCO2= expData(expIndex).listIDsInCueCO2(IDsToKeep);
        expData(expIndex).listTimeInCueCO2= expData(expIndex).listTimeInCueCO2(IDsToKeep);
        expData(expIndex).listTrajTimeCO2= tmpListTrajTime(IDsToKeep);
    end
end
clear tmpIndex tmpListTrajTime
%toc;


%mosqType=['m0';'m1';'wt';'l1';'l2';'l3';'l4';'l5';'l6'];
mosqTypes={'m0';'m0';'m1';'wt';'l1';'l4';'l3';'l5';'l2';'l6';'l1';'l6';'l2';'l5';'wt';'l3';'l4';'l5';'l3';'l4';'wt'};
uniqueTypes={'m0';'m1';'wt';'l1';'l2';'l3';'l4';'l5';'l6'};
figure()
for typeIndex=1:size(uniqueTypes,1)
    % Group experiments data in function of type of mosquitoes
    %typeValue= uniqueTypes(typeIndex);
    %i= find(strcmp(mosqTypes, typeValue));
    i= [3;4;5;6];
    timeTrajNearCueAIR=[];
    timeTrajNearCueCO2=[];
    for expIndex=i'
        % FOR AIR
        tmpNorm= bsxfun(@rdivide, (expData(expIndex).listTimeInCueAIR*100),expData(expIndex).listTrajTimeAIR);
        timeTrajNearCueAIR= [timeTrajNearCueAIR, tmpNorm];
        % Add percentage values to the strucutre
        expData(expIndex).percentNearCueAIR= tmpNorm;
        %FOR CO2
        tmpNorm= bsxfun(@rdivide, (expData(expIndex).listTimeInCueCO2*100),expData(expIndex).listTrajTimeCO2);
        timeTrajNearCueCO2= [timeTrajNearCueCO2, tmpNorm];
        % Add percentage values to the strucutre
        expData(expIndex).percentNearCueCO2= tmpNorm;
    end
    
    %Save values in table
    timeInCuePerTypeNorm(typeIndex).type= typeValue; 
    timeInCuePerTypeNorm(typeIndex).timeTrajNearCueAIR= timeTrajNearCueAIR;
    timeInCuePerTypeNorm(typeIndex).timeTrajNearCueCO2= timeTrajNearCueCO2;
    %Generate a SCATTER PLOT for this type of mutation
    subplot(3,3,typeIndex)
    scatter(ones(size(timeTrajNearCueAIR,2),1), timeTrajNearCueAIR, 10, 'jitter', 'on', 'MarkerFaceColor', [0, 0.5, 0], 'MarkerFaceAlpha', 0.5);
    hold on
    scatter(ones(size(timeTrajNearCueCO2,2),1)*2, timeTrajNearCueCO2, 10, 'jitter', 'on', 'MarkerFaceColor', [0.5, 0, 0], 'MarkerFaceAlpha', 0.5);
    % plot also the means values 
    y= [mean(timeTrajNearCueAIR), mean(timeTrajNearCueAIR)];
    plot([0.9, 1.1], y, 'g', 'LineWidth',5,'MarkerEdgeColor',[0 0.5 0]);
    y= [mean(timeTrajNearCueCO2), mean(timeTrajNearCueCO2)];
    plot([1.9, 2.1], y, 'r', 'LineWidth',5,'MarkerEdgeColor',[0.5 0 0]);
    hold off
    lg= [{' With AIR'}, {'With CO2'}];
    %lg= [{' With AIR'}, {'With CO2'}, {'With AIR (postCO2)'}];
    legend(lg);
    title(strcat('% time near black cue per trajectory (',typeValue,')'));
    xlim([0,3]);
    ylim([0, 100]);
    ylabel('% of time')
    %typeIndex= typeIndex+1;
end    

% ====================================
%=====================================





% ====================================
%=====================================

% ============= PrefIndex per 2x2cm2 square in XZ plane ===============
% From all the trajectories inside the cues, pick the ones we have X time before they entered the cue.
% Divide the XZ plane of the WT in 2x 2 cm2 squares (sub-square)
% “Generate the PI for each sub-squareas the amount of time spent on the side of the wind tunnel of the test object compared to the control object, divided by their sum.” (Each sub-square will have a list of PI values (one for each ID that pass around that cue)) and all the cubes crossed by a trajectory will have the same PI vale (for that trajectory ID)
% "Generate the mean PI value for each one of the sub-square accros all trajectories and its 95% confidence interval"

% Load for each trajectory ID its first time it enters the a volume
odorChecked='CO2';
if (exist('outputPath', 'var') &&  exist('outputFolder', 'var'))
    % Load files
    filesPath= strcat(outputPath, outputFolder, 'analysisData\');
    cd(filesPath);
    filesList=dir(strcat('*countsInsideCueVol_',odorChecked,'.xlsx'));
    cd(workspace);
    for fileIndex= 1:length(filesList)
        fileName= filesList(fileIndex).name;
        disp(strcat(' - Working with file: ', {' '}, fileName));
        filesName(fileIndex)= {fileName(1:15)};
        % data loaded structure: [vCuePosition	objID	timeStamp NaN initialTimeStamp finalTimeStamp]
        dataFromExcel = table2array(readtable(strcat(filesPath, fileName)));
        if isempty(dataFromExcel)
            disp(strcat('   -> No counts detected in file: ',{' '}, fileName));
            % The matrices p1..p4 and t1, t2 will keep its 0 values in row assigned to this empty file 
        else
            %Sort the data
            sortedData= sortrows(dataFromExcel(:,[1,2,3,6]));
            %find the first appearence of the position 2 
            split= find(sortedData(:,1) ==2,1);
            %Find the information related to the Pos1
            uniqueIDs= unique(sortedData(1:(split-1),2));
            %initialize the structure containing the data to use in this analysis
            expData(fileIndex).idListP1= uniqueIDs';
            expData(fileIndex).tsListP1= zeros(1,length(uniqueIDs));
            expData(fileIndex).firstIxP1= zeros(1,length(uniqueIDs));
            expData(fileIndex).lastIxP1= zeros(1,length(uniqueIDs));
            for i= 1:length(uniqueIDs)
                %Find the TimeStamp and Raw data index for the first time the ID enters the volume
                tempIndex= find(sortedData(1:(split-1),2)==uniqueIDs(i),1);
                expData(fileIndex).tsListP1(i)= sortedData(tempIndex, 3);
                expData(fileIndex).firstIxP1(i)= sortedData(tempIndex, 4);
                %Find the TimeStamp and Raw data index for the last time the ID is inside the volume
                tempIndex= find(sortedData(:,2)==uniqueIDs(i),1, 'last');
                expData(fileIndex).lastIxP1(i)= sortedData(tempIndex, 4);
            end 
            % Find the information for possition 2
            uniqueIDs= unique(sortedData(split:end,2));
            %initialize the structure containing the data to use in this analysis
            expData(fileIndex).idListP2= uniqueIDs';
            expData(fileIndex).tsListP2= zeros(1,length(uniqueIDs));
            expData(fileIndex).firstIxP2= zeros(1,length(uniqueIDs));
            expData(fileIndex).lastIxP2= zeros(1,length(uniqueIDs));
            for i= 1:length(uniqueIDs)
                % Find the TimeStamp and Raw data index for the first time the ID enters the volume and
                % extrapolate to the full sortedData set (instead of starting for the split position)
                tempIndex= find(sortedData(split:end,2)==uniqueIDs(i),1)+(split-1);
                expData(fileIndex).tsListP2(i)= sortedData(tempIndex, 3);
                expData(fileIndex).firstIxP2(i)= sortedData(tempIndex, 4);
                %Find the TimeStamp and Raw data index for the first time the ID enters the volume
                tempIndex= find(sortedData(:,2)==uniqueIDs(i),1, 'last');
                expData(fileIndex).lastIxP2(i)= sortedData(tempIndex, 4);
            end 
        end
    end
end
clear dataFromExcel dataSorted uniqueIDs fstIndex

% The expData structure has been created witht he following fields (each row== experiment date):
% - idList: List of IDs detected inside one of the volume
% - tsList: list with the first time the ID was inside the volume (IDs in idList)
% - firstIxList: index related to the RAW DATA that can be loaded in the
%           "dataset" structure for the moment the ID was inside a volume
% - lasttIxList: index related to the RAW DATA that can be loaded in the
%           "dataset" structure for the moment the ID left a volume

% FOR TEST
%  expData2= expData;
%  expData(1).idListP1= [90,94,110,111,112,113,115,116,117,120,130];
%  expData(1).tsListP1= [1,10,20,29,43,32,10,40,50,55, 60];
%  expData(1).idListP2= [90,99,105,110,115,117,120,123,130];
%  expData(1).tsListP2= [5,20,21,10,30,39,57,59,59];

%Find if any ID has visit both cues, then keep only the cue visited as first choice
%expIndex=1;
for expIndex=1:length(expData)
    disp(strcat('# of Ids visiting P1: ', num2str(length(expData(expIndex).idListP1))));
    disp(strcat('# of Ids visiting P2: ', num2str(length(expData(expIndex).idListP2))));
    visitBothCues= ismember(expData(expIndex).idListP2,expData(expIndex).idListP1);
    % visitBothCues==0 ID not repeated and visitBothCues==1 ID repeated in both cues
    if any(visitBothCues)
        [repeats,value]= groupcounts(visitBothCues');    %how many 0 and 1 are repeated
        disp(strcat('# of Ids visiting both cues: ', num2str(repeats(2))));
        for i= find(visitBothCues ==1)
            %Check which of both cues was visited as first choice
            j= find(expData(expIndex).idListP1 == expData(expIndex).idListP2(i));
            %j= find(ctrIdP1 == ctrIdP2(i));
            %if (ctrTsP1(j) - ctrTsP2(i)) > 0
            if (expData(expIndex).tsListP1(j) -  expData(expIndex).tsListP2(i)) > 0
                %The visit in Pos1 happened after the visit of Pos2
                % Change the information for ID to 0 to remove later from matrix
                expData(expIndex).idListP1(j)= 0;
                expData(expIndex).tsListP1(j)= 0;
                expData(expIndex).firstIxP1(j)= 0;
                expData(expIndex).lastIxP1(j)= 0;
            else
                %The visit in Pos2 happened after the visit of Pos1             
                % Change the information for ID to 0 to remove later from matrix
                expData(expIndex).idListP2(i)= 0;
                expData(expIndex).tsListP2(i)= 0;
                expData(expIndex).firstIxP2(i)= 0;
                expData(expIndex).lastIxP2(i)= 0;
            end
        end
    end
end
%Check how many value will be eleminated from each idList
disp(strcat('# of Ids visiting P1 as second cue choice: ', num2str(numel(find(expData(1).idListP1==0)))));
disp(strcat('# of Ids visiting P2 as second cue choice: ', num2str(numel(find(expData(1).idListP2==0)))));

%Clean all values for the IDs set to 0
expData(expIndex).idListP1(expData(expIndex).idListP1 == 0) =[];
expData(expIndex).tsListP1(expData(expIndex).tsListP1 == 0) =[];
expData(expIndex).firstIxP1(expData(expIndex).firstIxP1 == 0)= [];
expData(expIndex).lastIxP1(expData(expIndex).lastIxP1 == 0)= [];
expData(expIndex).idListP2(expData(expIndex).idListP2 == 0) =[];
expData(expIndex).tsListP2(expData(expIndex).tsListP2 == 0) =[];
expData(expIndex).firstIxP2(expData(expIndex).firstIxP2 == 0)= [];
expData(expIndex).lastIxP2(expData(expIndex).lastIxP2 == 0)= [];

clear repeat values visitBothCues i j expIndex


% Generate the sqSz x sqSz cm2 square map of the xz plane
% sqSz= 0.02; 
% sqInZ= round(lim_z/sqSz);
% % Only half X axis
% %sqInX= round(lim_x/sqSz);
% %Full X axis
% sqInX= round(2*lim_x/sqSz);
% sqInY= round(2*lim_y/sqSz);


% Now we will need to work with the RAW DATA to look for the IDs that we
% have recored longer than 2 seconds

% %Check RAW DATA indexes (.ixList) is linked to the same IDs as .idList
% for expIndex= 1: length(expData)
%     tempIDs= dataset(expIndex).attr_id(expData(expIndex).firstIxP1);
%     kk= expData(expIndex).idListP1 - tempIDs';
%     if sum(kk) ~= 0
%         disp(strcat('There is an error with data from:', filesName(expIndex)));
%     end
% end



% % ========== Spatial representation of preference index prior vue visit ===============
% Look for the duration in time for each trajectory ID and keep all the
% trajecotories with tLimit sec or more before the trajecotry enters a volume
tLimit= 2; % num opf seconds before the trajectory reach teh volume
mainIndex= 1;
%Initialize parameter for multidimensional array
xEdges= -lim_x:0.02:lim_x; %-0.76:0.02:0;
zEdges= 0:0.02:lim_z;
XZplane= NaN(length(xEdges)-1, length(zEdges)-1, 1);

for expIndex=1:length(expData)
    count= 0;
    for i= 1:length(expData(expIndex).idListP1)
        %Find the indexes for the choosen trajectory ID
        ix= find(dataset(expIndex).attr_id == expData(expIndex).idListP1(i));
        if (expData(expIndex).tsListP1(i) - dataset(expIndex).attr_time(ix(1))) >= tLimit
            %Find the first point for the trajectory detected 2 sec before the traj reach the volume
            trajInitialTime= expData(expIndex).tsListP1(i)-tLimit;
            trajInitialIndex= find(dataset(expIndex).attr_time(ix) >= trajInitialTime, 1);
            % extrapolate the index from the ix index List to the real position omn the RAW DATA
            indexesToUse= ix(trajInitialIndex):expData(expIndex).lastIxP1(i);

            % Find how many points from trajectory are in each of the sides of the Y axis
            % and generate the PI value to be added for each square in the XZ plane
            ctsP1= numel(find(dataset(expIndex).attr_y(indexesToUse) > 0));
            ctsP2= numel(find(dataset(expIndex).attr_y(indexesToUse) < 0));
            if baseColorIndexList(expIndex) == 1
                %Base cue was placed in the Pos-1
                tempPI= (ctsP2-ctsP1)/(ctsP1+ctsP2);
            else
                %Base cue was placed in the Pos-2
                tempPI= (ctsP1-ctsP2)/(ctsP1+ctsP2);
            end

            %Create the XZ plane as a square grid and find which squares
            %where visited by trajectory
            h=histcounts2(dataset(expIndex).attr_x(indexesToUse),dataset(expIndex).attr_z(indexesToUse),xEdges, zEdges);
            
            %Find in which bins were visited by the trajectory ID
            p= h > 0;
            p2= h== 0;
            %Store the PI value generated for this trajectory ID in the bins visited by it
            h(p)= tempPI;
            h(p2)=NaN;
            %Add the histcounts bins plane to the multidimensional array (3rd dimmension) as a new page
            XZplane(:,:,mainIndex)= h;
            mainIndex= mainIndex+1;
            count= count+1;
        end
    end
    %disp(strcat('total trajectories visitn Pos-1 as first choice longer that 2 sec', {' '}, num2str(count)));
    count= 0;
    for i= 1:length(expData(expIndex).idListP2)
        ix= find(dataset(expIndex).attr_id == expData(expIndex).idListP2(i));
        if (expData(expIndex).tsListP2(i) - dataset(expIndex).attr_time(ix(1))) >= tLimit
            %Find the first point for the trajectory detected 2 sec before the traj reach the volume
            trajInitialTime= expData(expIndex).tsListP2(i)-tLimit;
            trajInitialIndex= find(dataset(expIndex).attr_time(ix) >= trajInitialTime, 1);
            % extrapolate the index from the ix index List to the real position omn the RAW DATA
            indexesToUse= ix(trajInitialIndex):expData(expIndex).lastIxP2(i);

            % Find how many points from trajectory are in each of the sides of the Y axis
            % and generate the PI value to be added for each square in the XZ plane
            ctsP1= numel(find(dataset(expIndex).attr_y(indexesToUse) > 0));
            ctsP2= numel(find(dataset(expIndex).attr_y(indexesToUse) < 0));
            if baseColorIndexList(expIndex) == 1
                %Base cue was placed in the Pos-1
                tempPI= (ctsP2-ctsP1)/(ctsP1+ctsP2);
            else
                %Base cue was placed in the Pos-2
                tempPI= (ctsP1-ctsP2)/(ctsP1+ctsP2);
            end

            %Create the XZ plane as a square grid and find which squares where visited by trajectory
            h=histcounts2(dataset(expIndex).attr_x(indexesToUse),dataset(expIndex).attr_z(indexesToUse),xEdges, zEdges);
            
            %Find in which bins were visited by the trajectory ID
            p= h > 0;
            p2= h==0;
            %Store the PI value generated for this trajectory ID in the bins visited by it
            h(p)= tempPI;
            h(p2)=NaN;
            % Add the trajectory ID information to the multidimensional
            % array (3rd dimension is ID index)
            %Add the histcounts bins plane to the multidimensional array (3rd dimmension)           
            XZplane(:,:,mainIndex)= h;
            mainIndex= mainIndex + 1;
            count= count + 1;
        end
    end
    disp(strcat('total trajectories visitn Pos-2 as first choice longer that 2 sec',{' '}, num2str(count)));  
end
clear expIndex tLimit tempPI h p count ctsP1 ctsP2 indexesToUse trajInitialtime trajInitialIndex ix i p1Counts

% Generate the mean PI value for each bin in the XZplane (use 3rd dimension)
meanXZplane= nanmean(XZplane,3);
% plot the mean PI values for eacjh bin in the XZ plane
figure();
minv = min(min(meanXZplane));
maxv = max(max(meanXZplane));
meanXZplane(isnan(meanXZplane)) = minv-((maxv-minv)/max(size(meanXZplane)));
ddd=[0 0 0;hot(10)];
colormap(ddd);
pcolor(meanXZplane')
cb= colorbar;
set(cb, 'ylim', [-1 1])
xlabel('X Axis');
ylabel('Y Axis');
title('Spatial representation of preference index prior vue visit');
clear cb


% figure();
% pcolor(meanXZplane','ShowEmptyBins','on')
% colormap(hot)
% cb= colorbar;
% set(cb, 'ylim', [-1 1])
% xlabel('X Axis');
% ylabel('Y Axis');
% title('Spatial representation of preference index prior vue visit');
% clear cb
%==================================================
%==================================================





% 
% %Find if any ID has visit both cues, then keep only the cue visited as
% %first choice
% expIndex=1;
% for expIndex=1:length(expData)
%     count= 0;
%     countP1=0;
%     countP2=0;
%     tempIndexesI=[];
%     disp(strcat('# of Ids visiting P1: ', num2str(length(expData(expIndex).idListP1))));
%     disp(strcat('# of Ids visiting P2: ', num2str(length(expData(expIndex).idListP2))));
%     visitBothCues= ismember(expData(expIndex).idListP2,expData(expIndex).idListP1);
%     % visitBothCues==0 ID not repeated and visitBothCues==1 ID repeated in both cues
%     if any(visitBothCues)
%         actrIdP1= expData(expIndex).idListP1;
%         ctrTsP1= expData(expIndex).tsListP1;
%         ctrIdP2= expData(expIndex).idListP2;
%         ctrTsP2= expData(expIndex).tsListP2;
%         ctrIxfP1= expData(expIndex).firstIxP1;
%         ctrIxlP1= expData(expIndex).lastIxP1;
%         ctrIxfP2= expData(expIndex).firstIxP2;
%         ctrIxlP2= expData(expIndex).lastIxP2;
%         [repeats,value]= groupcounts(visitBothCues');    %how many 0 and 1 are repeated
%         disp(strcat('# of Ids visiting both cues: ', num2str(repeats(2))));
%         for i= find(visitBothCues ==1)
%             tempIndexesI= [tempIndexesI, i]; 
%             count=count+1;
%             %Check which of both cues was visited as first choice
%             %j= find(expData(expIndex).idListP1 == expData(expIndex).idListP2(i));
%             j= find(ctrIdP1 == ctrIdP2(i));
%             if (ctrTsP1(j) - ctrTsP2(i)) > 0
%                 countP1=countP1+ 1;
%                 %The visit in Pos1 happened after the visit of Pos2
%                 if (j==1)
%                     % Delete the first ID and its ts and indexes from the expData set
%                     expData(expIndex).idListP1= expData(expIndex).idListP1((j+1):end);
%                     expData(expIndex).tsListP1= expData(expIndex).tsListP1((j+1):end);
%                     expData(expIndex).firstIxP1= expData(expIndex).firstIxP1((j+1):end);
%                     expData(expIndex).lastIxP1= expData(expIndex).lastIxP1((j+1):end);
%                 elseif (j > 1) && (j < length(ctrIdP1))
%                     % Delete the ID and its ts and indexes from the expData set
%                     expData(expIndex).idListP1= [expData(expIndex).idListP1(1:(j-1)),  ctrIdP1((j+1):end)];
%                     expData(expIndex).tsListP1= [expData(expIndex).tsListP1(1:(j-1)),  ctrTsP1((j+1):end)];
%                     expData(expIndex).firstIxP1=[expData(expIndex).firstIxP1(1:(j-1)), ctrIxfP1((j+1):end)];
%                     expData(expIndex).lastIxP1= [expData(expIndex).lastIxP1(1:(j-1)),  ctrIxlP1((j+1):end)];
%                 elseif(j == length(ctrIdP1))
%                     %Delete form the expData the information for last ID
%                     expData(expIndex).idListP1= expData(expIndex).idListP1(1:(end-1));
%                     expData(expIndex).tsListP1= expData(expIndex).tsListP1(1:(end-1));
%                     expData(expIndex).firstIxP1= expData(expIndex).firstIxP1(1:(end-1));
%                     expData(expIndex).lastIxP1= expData(expIndex).lastIxP1(1:(end-1));
%                 end
%             else
%                 countP2= countP2+1;
%                 %The visit in Pos2 happened after the visit of Pos1             
%                 if (i==1)
%                     % Delete the first ID and its ts and indexes from the expData set
%                     expData(expIndex).idListP2= expData(expIndex).idListP2((i+1):end);
%                     expData(expIndex).tsdListP2= expData(expIndex).tsListP2((i+1):end);
%                     expData(expIndex).firstIxP2= expData(expIndex).firstIxP2((i+1):end);
%                     expData(expIndex).lastIxP2= expData(expIndex).lastIxP2((i+1):end);
%                 elseif ((i > 1) && (i < length(ctrIdP2)))
%                     % Delete the ID and its ts and indexes from the expData set
%                     expData(expIndex).idListP2= [expData(expIndex).idListP2(1:(i-1)),  ctrIdP2((i+1):end)];
%                     expData(expIndex).tsListP2= [expData(expIndex).tsListP2(1:(i-1)),  ctrTsP2((i+1):end)];
%                     expData(expIndex).firstIxP2=[expData(expIndex).firstIxP2(1:(i-1)), ctrIxfP2((i+1):end)];
%                     expData(expIndex).lastIxP2= [expData(expIndex).lastIxP2(1:(i-1)),  ctrIxlP2((i+1):end)];
%                 elseif(i == length(ctrIdP2))
%                     %Delete form the expData the information for last ID
%                     expData(expIndex).idListP2= expData(expIndex).idListP2(1:(i-1));
%                     expData(expIndex).tsListP2= expData(expIndex).tsListP2(1:(i-1));
%                     expData(expIndex).firstIxP2= expData(expIndex).firstIxP2(1:(i-1));
%                     expData(expIndex).lastIxP2= expData(expIndex).lastIxP2(1:(i-1));
%                 end 
%             end
%         end
%         count= 0;
%     end
% end
% disp(strcat('# of Ids visiting P1, firstChoice only: ', num2str(length(expData(expIndex).idListP1))));
% disp(strcat('# of Ids visiting P2, firstChoice only: ', num2str(length(expData(expIndex).idListP2))));
% disp(count)
% clear ctrIdP1 ctrTsP1 ctrIxfP1 ctrIxlP1
% clear ctrIdP2 ctrTsP2 ctrIxfP2 ctrIxlP2







% ============================
% ============================



% ============= PrefIndex per trajectory ID ===============
% Initialize the matrices for the timestamps
% Set the amount of groups to divide the CO2 part experiment
numGrps=120; %=120 =60 sec per group || =4 =1800 sec per grp
for fileIndex= 1:length(filesList)
    fileName= filesList(fileIndex).name;
    %xNames(fileIndex)= cellstr(filesList(fileIndex).name(1:8));
    disp(strcat(' - Working with file: ', {' '}, fileName));
    filesName(fileIndex)= {fileName(1:8)};
    % data loaded structure: [vCuePosition	objID	timeStamp NaN initialTimeStamp finalTimeStamp]
    dataFromExcel = table2array(readtable(strcat(filesPath, fileName)));
    if isempty(dataFromExcel)
        disp(strcat('   -> No counts detected in file: ',{' '}, fileName));
        % The matrices p1..p4 and t1, t2 will keep its 0 values in row assigned to this empty file 
    else
        initialTS= dataFromExcel(1,4);    
        finalTS= dataFromExcel(1,5);
        grpSz= (finalTS - initialTS)/numGrps;
        grpThreshold= zeros(1, numGrps+1);
        for i=0:numGrps
            grpThreshold(i+1)= initialTS+(grpSz*i);
        end       
        sortedData= sortrows(dataFromExcel(:,1:3));
        %find the first appearence of the position 2 
        split= find(sortedData(:,1) ==2,1);

        % ** Counts the number of repetiions for each ID in each cue **
        % - Cue in Pos 1
        [IDsCtsP1,IDsNameP1] = groupcounts(sortedData(1:(split-1),2));
        % - Cue in Pos 2
        [IDsCtsP2,IDsNameP2] = groupcounts(sortedData(split:end,2));
        i=1;
        uniqueID= unique(sortedData(:,2));
        for id=uniqueID'
            iP1= find(IDsNameP1 == id);
            iP2= find(IDsNameP2 == id);
            if any(iP1) && any(iP2)
                % ID has been near both cues
                if baseColorIndexList(fileIndex) == 1
                    tempPI= (IDsCtsP2(iP2) - IDsCtsP1(iP1)) / (IDsCtsP2(iP2) + IDsCtsP1(iP1));
                else
                    tempPI= (IDsCtsP1(iP1) - IDsCtsP2(iP2)) / (IDsCtsP2(iP2) + IDsCtsP1(iP1));    
                end
            elseif isempty(iP2)
                % ID has been near POS 1 only
                if baseColorIndexList(fileIndex) == 1
                    tempPI= (0 - IDsCtsP1(iP1)) / IDsCtsP1(iP1);
                else
                    tempPI= (IDsCtsP1(iP1) - 0) / IDsCtsP1(iP1);    
                end
            elseif isempty(iP1)
                % ID has been near POS 2 only
                if baseColorIndexList(fileIndex) == 1
                    tempPI= (IDsCtsP2(iP2) - 0) / IDsCtsP2(iP2);
                else
                    tempPI= (0 - IDsCtsP2(iP2)) / IDsCtsP2(iP2);    
                end
            end
            piPerID_list(i)= tempPI;
            i= i+1;
        end
        piPerID_list(isnan(piPerID_list))= [];
        piCounts_list(fileIndex).mean= mean(piPerID_list);
        piCounts_list(fileIndex).piPerID= piPerID_list;
        
        % ** PI for ID looking at the time spent by insect near each cue **
        i=1;
        uniqueID= unique(sortedData(:,2));
        for id=uniqueID'
            iP1= find(sortedData(1:(split-1),2) == id);
            iP2= find(sortedData(split:end,2) == id);
            if any(iP1)
                % ID has been near POS 1
                % calculate the Delta time between timestamps values 
                tsDiffP1= diff(sortedData(iP1, 3));
                % find which  sequential counts didn't happen "~=consecutively"
                outList= find(tsDiffP1 > 0.05);
                t=0;
                k=1;
                if any(outList)
                    for outValue= outList'
                        % the insect left the volume at a given moment
                        t= t + sum(tsDiffP1(k:(outValue-1)));
                        k= outValue+1;
                    end
                    t= t + sum(tsDiffP1(k:end));
                    % Another posibility is to add the minimum time for the isolated detectections (outList)
                    %Add the minimum time (given by the fps) for the isolated detectections (outList)
                    %t= t + length(outList)*1/fps;                 
                else
                    %the insects has been inside the volume all the time
                    t= sum(tsDiffP1);
                end
                timeInP1(i)= t;
            else
                % the ID has not been near POS 2                
                timeInP1(i)= 0;
            end
            if any(iP2)
                % ID has been near POS 2
                tsDiffP2= diff(sortedData(iP2+(split-1), 3));
                % find which  sequential counts didn't happen "~=consecutively"
                outList= find(tsDiffP2 > 0.05);
                t=0;
                k=1;
                if any(outList)
                    for outValue= outList'
                        % the insect left the volume at a given moment
                        t= t + sum(tsDiffP2(k:(outValue-1)));
                        k= outValue+1;
                    end
                    t= t + sum(tsDiffP2(k:end));
                    % Another posibility is to add the minimum time for the isolated detectections (outList)                    
                    %Add the minimum time (given by the fps) for the isolated detectections (outList)
                    %t= t + length(outList)*1/fps;                        
                else
                    %the insects has been inside the volume all the time
                    t= sum(tsDiffP2);            
                end
                timeInP2(i)= t;
            else
                % the ID has not been near POS 2
                timeInP2(i)=0;
            end
            
            %Generate the PI for the given ID
            if baseColorIndexList(fileIndex) == 1
                tempPI= (timeInP2(i) - timeInP1(i)) / (timeInP1(i) + timeInP2(i));
            else
                tempPI= (timeInP1(i) - timeInP2(i)) / (timeInP1(i) + timeInP2(i));    
            end
            piPerID_list(i)= tempPI;
            i= i+1;
        end      
        piPerID_list(isnan(piPerID_list))= [];      
        piTime_list(fileIndex).mean= mean(piPerID_list);
        piTime_list(fileIndex).piPerID= piPerID_list;
    end
end

% Ploting Pref Index per counts for each ID near the cues
figure()
for k=1:length(filesList)
    subplot(2,length(filesList),k)
    scatter(ones(1, size(piCounts_list(k).piPerID,2)),piCounts_list(k).piPerID, 'jitter', 'on')
    hold on
    plot(median(piCounts_list(k).piPerID), 'ro', 'MarkerSize',10,'MarkerEdgeColor','r', 'MarkerFaceColor',[0.9, 0.5, 0.5]);
    y= [mean(piCounts_list(k).piPerID), mean(piCounts_list(k).piPerID)];
    plot([0.7, 1.3], y, '--k', 'MarkerSize',15,'MarkerEdgeColor','k');
    hold off
    title('PI for counts per ID');
    ylabel('Preference Index'); 
    xlabel(strcat('# of IDs= ',num2str(length(piCounts_list(k).piPerID))));
    subplot(2,5,k+length(filesList))
    histogram(piCounts_list(k).piPerID, [-1:0.1:1]);
    ylabel(' Repeats for value');
    xlabel('Preference Index')
    title(filesName(k))
end
    
% Ploting Pref Index per time spent each ID near the cues
figure()
for k=1:length(filesList)
    subplot(2,length(filesList),k)
    scatter(ones(1, size(piTime_list(k).piPerID,2)),piTime_list(k).piPerID, 'jitter', 'on')
    hold on
    plot(median(piTime_list(k).piPerID), 'ro', 'MarkerSize',10,'MarkerEdgeColor','r', 'MarkerFaceColor',[0.9,0.5,0.5]);
    y= [mean(piTime_list(k).piPerID), mean(piTime_list(k).piPerID)];
    plot([0.7, 1.3], y, '--k', 'MarkerSize',15,'MarkerEdgeColor','k');
    hold off
    title('PI for time spent in cue per ID');
    ylabel('Preference Index'); 
    xlim([0.5, 1.5]);
    xlabel(strcat('# of IDs= ',num2str(length(piTime_list(k).piPerID))));
    subplot(2,length(filesList),k+length(filesList))
    histogram(piTime_list(k).piPerID, [-1:0.1:1]);
    ylabel(' Repeats for value');
    xlabel('Preference Index')
    title(filesName(k))
end
% ==================================================
% =================================================    

% ======= Comprare scatter plot with time spent per time group of 60 sec
duration= 60;
totalExp= length(filesNameList);
% Find in which position the base and test color wehre placed
p1BaseC= find(baseColorIndexList == 1);
p2BaseC= find(baseColorIndexList== 2);
sumTimeInBaseClrCO2= cumsum([t1CO2(p1BaseC, 1:duration); t2CO2(p2BaseC, 1:duration)], 2);
sumTimeInTestClrCO2= cumsum([t2CO2(p1BaseC, 1:duration); t1CO2(p2BaseC, 1:duration)], 2);
kkBC1= [];
kkTC1= [];
kkBC2= [];
kkTC2= [];
testInBC=[];
testInTC=[];
for fileIndex= 1:totalExp
    %Group the data related to the Base and Tested cues
    if baseColorIndexList(fileIndex) == 1    
        kkBC1= [kkBC1; t1CO2(fileIndex, :)];
        kkTC2= [kkTC2; t2CO2(fileIndex, :)];
        testInBC= [testInBC; t1CO2(fileIndex, :)];
        testInTC= [testInTC; t2CO2(fileIndex, :)];
    else
        kkBC2= [kkBC2; t2CO2(fileIndex, :)];
        kkTC1= [kkTC1; t1CO2(fileIndex, :)];
        testInBC= [testInBC; t2CO2(fileIndex, :)];
        testInTC= [testInTC; t1CO2(fileIndex, :)];
    end
end
% generate cumsums and plot with cumsums above
sumkkBC1= cumsum(kkBC1, 2);
sumkkTC1= cumsum(kkTC1, 2);
sumkkBC2= cumsum(kkBC2, 2);
sumkkTC2= cumsum(kkTC2, 2);
sumtestInBC=cumsum(testInBC, 2);
sumtestInTC=cumsum(testInTC, 2);

figure()
subplot(2,1,1)
plot(mean([sumkkBC1; sumkkBC2]))
hold on
plot(mean(sumtestInBC))
plot(mean(sumTimeInBaseClrCO2), 'r');
hold off
ylim([0, max(sumTimeInTestClrCO2(:,end))])
lg=[{'kkBC'}, {'testInBC'}, {'sumTimeinBaseClrCO2'}];
legend(lg);
subplot(2,1,2)
plot(mean([sumkkTC2; sumkkTC1])')  
hold on
plot(mean(sumtestInTC)')
plot(mean(sumTimeInTestClrCO2'));
hold off
ylim([0, max(sumTimeInTestClrCO2(:,end))])
lg=[{'kkTC'}, {'testInTC'}, {'sumTimeInTestClrCO2'}];
legend(lg);

    
    
    

    
    
        
        
        
        
        
        

for grpIndex= 1:numGrps %length(grpThreshold)
    %Select the indexes from the data that verify the time group
    %conditions
    dataIndexes= find(dataFromExcel(:,3) < grpThreshold(grpIndex+1) & dataFromExcel(:,3) >= grpThreshold(grpIndex));

    p1(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 1);
    p2(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 2);
    p3(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 3);
    p4(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 4);
    % To load timestamp data (time spent by insect near thew cue
    % Sorted the subset of data by cue position
    sortedData= sortrows(dataFromExcel(dataIndexes,1:3));
    %find the first appearence of the position 2 
    split= find(sortedData(:,1) ==2,1);

    % calculate the Delta time between timestamps values 
    tsDiffP1= diff(sortedData(1:(split-1), 3));
    tsDiffP2= diff(sortedData(split:end, 3));
    % For POS-1
    % find which  sequential counts didn't happen "~=consecutively"
    out= find(tsDiffP1 > 0.05);
    t=0;
    k=1;
    if any(out)
        for i= out
            % the insect left the volume at a given moment
            t= t + sum(tsDiffP1(k:(i-1)));
            k= i+1;
        end
        t= t + sum(tsDiffP1(k:end));
    else
        %the insects has been inside the volume all the time
        t= sum(tsDiffP1);
    end
    t1(fileIndex, grpIndex)= t;
    % For POS-2
    % find which  sequential counts didn't happen "~=consecutively"
    out= find(tsDiffP2 > 0.05);
    t=0;
    k=1;
    if any(out)
        for i= out
            % the insect left the volume at a given moment
            t= t + sum(tsDiffP2(k:(i-1)));
            k= i+1;
        end
        t= t + sum(tsDiffP2(k:end));
    else
        %the insects has been inside the volume all the time
        t= sum(tsDiffP2);
    end
    t2(fileIndex, grpIndex)= t;

end

















%==========================================================


%(amount of time trajs are in cue_1 - amount of time trajs are in cue_2) /(total amount of times in a trajectory)
%Test dataset
% a= [1 ,10; 1, 11; 1 ,12; 1, 13; 1 ,14; 1, 15;1 ,16; 1, 17; ...
%     2 ,11; 2, 13; 2 ,15; 2, 17; 2 ,19; 3, 15;1 ,26; 1, 27; ...
%     3 ,21; 3, 23; 3 ,25; 2, 21; 2 ,29; 3, 45;4 ,16; 4, 17];
% cuesP= randi(2,length(a(:,1)),1);
% a= horzcat(cuesP, a)

%a.txt contains a matrix to use in test


%sort matrix in fct of the cue 
sortedData= sortrows(a);
%find the first appearence of the position 2 
split= find(sortedData(:,1) ==2,1);
%counts the number of repetiions for each ID in each cue
%Cue in Pos 1
[IDsCtsP1,IDsGrpP1] = groupcounts(sortedData(1:(split-1),2));
%Cue in Pos 2
[IDsCtsP2,IDsGrpP2] = groupcounts(sortedData(split:end,2));
%Find the estimation of the mean duration any insect spent in each cue
durInPos1= mean([IDsCtsP1; IDsCtsP2])/evalin('base',fps);
    


% TIMESTAMP INSTANT NO LOADED FROM XLS FILE    
    filesPath= strcat(outputPath, outputFolder);
    numGrps=120; %=120 =60 sec per group || =4 =1800 sec per grp
    fileIndex= 3;
    fileName= strcat(dataset(fileIndex).fileName(1:15),'_countsInsideCueVol_CO2.xlsx');
    
    %Matrices with the different ts values [(t1-t0), ..., (tn- t(n-1))] 
    tDiff1 = zeros(length(filesList),numGrps); 
    tDiff2 = zeros(length(filesList),numGrps);
       


    fileName= filesList(fileIndex).name;
    disp(strcat(' - Working with file: ', {' '}, fileName));
    %filesName(fileIndex)= {fileName(1:8)};
    % data loaded structure: [vCuePosition	objID	timeStamp NaN initialTimeStamp finalTimeStamp]
    dataFromExcel = table2array(readtable(strcat(filesPath, fileName)));
    initialTS= dataFromExcel(1,4);    
    finalTS= dataFromExcel(1,5);

    %sort matrix in fct of the cue-id-timestamp
    sortedData= sortrows(dataFromExcel(:,1:3));
    %find the first appearence of the position 2 
    split= find(sortedData(:,1) ==2,1);

    %Analyze the difference in between timestamps values (total, pos1, pos2)
    % This can show changes between Cues, IDs, or an ID leavin a
    % volumen and then coming back inside the same volume
    tsDiff= diff(sortedData(:,3))
    bar(tsDiff);
    ylabel('Delta values ')
    xlabel('Timestamps index');
    title(' Delta between timestamps of consecutive counts');
    figure()
    subplot(2,1,1)
    bar(tsDiff(1:split-2))
    ylabel('Delta values ')
    xlabel('Timestamps index');
    title(' Delta between timestamps of consecutive counts (Pos-1)');
    subplot(2,1,2)
    bar(tsDiff(split-1:end))
    ylabel('Delta values ')
    xlabel('Timestamps index');
    title(' Delta between timestamps of consecutive counts (Pos-1)');
    %mean diffrence between timestamps
    zm=mean(tsDiff);
    %Total number of differences between timestamps THAT FOLLOWS 0<tsDiff<0.3
    zm1=sum((tsDiff >0 & tsDiff <0.05));
    % '%' of differences smaller than 0.3 s
    zm2=(zm1*100)/length(tsDiff);
    disp(strcat(' * mean value for tsDiff (difference between timestamps):',num2str(zm)));
    disp(strcat(' * total # of tsDiff values 0>tsDiff(i)<0.3 :',num2str(zm1), '% of :',num2str(zm2)));
    zm1= sum((abs(tsDiff) >0 & abs(tsDiff) <0.3));
    zm2=(zm1*100)/length(tsDiff);
    disp(strcat(' * total # of tsDiff values 0>|tsDiff(i)|<0.3 :',num2str(zm1), '% of :',num2str(zm2)));

    %counts the number of repetiions for each ID in each cue
    %Cue in Pos 1
    [IDsCtsP1,IDsGrpP1] = groupcounts(sortedData(1:(split-1),2));
    %Cue in Pos 2
    [IDsCtsP2,IDsGrpP2] = groupcounts(sortedData(split:end,2));
    %Find the estimation of the total duration spent by all insects in each cue
    timeInPos1= sum(IDsCtsP1)/fps
    timeInPos2= sum(IDsCtsP2)/fps 
    % preference index towards POS 1
    prefIndexPos=(timeInPos1 - timeInPos2)/ (timeInPos1+timeInPos2)

    

    
        
        % === ESTIMATES HOW MUCH TIME AN INSECT SPEND INSIDE THE VOLUME
        %How much REAL time spend the insect selected in IDsToLook in the volume
        % In pos 1= idToLook(1)
        indexID= find(sortedData(:,2) == idToLook(1));
        %Check between first and last value
        durationEst= sortedData(indexID(end),3)- sortedData(indexID(1),3);     
        %Check its delta over ts
        auxDiff=diff(sortedData(indexID,3));
        out= find(auxDiff > 0.05);
        t=0;
        k=1;
        if any(out)
            for i= out
                % the insect left the volume at a given moment
                t= t + sum(auxDiff(k:(i-1)));
                k= i+1;
            end
            t= t + sum(auxDiff(k:end));
        else
            t= sum(auxDiff);
        end
        disp(strcat(' * 1st time estimation inside the volume: ',num2str(durationEst)));
        disp(strcat(' * 2nd time estimation inside the volume (ts delta <= 0.3): ',num2str(t)));
        
        % In pos 2= idToLook(2)
        indexID= find(sortedData(:,2) == idToLook(2));
        %Check between first and last value
        durationEst= sortedData(indexID(end),3)- sortedData(indexID(1),3);     
        %Check its delta over ts
        auxDiff=diff(sortedData(indexID,3));
        out= find(auxDiff > 0.3);
        t=0;
        k=1;
        if any(out)
            for i= out
                % the insect left the volume at a given moment
                t= t + sum(auxDiff(k:(i-1)));
                k= i+1;
            end
            t= t + sum(auxDiff(k:end));
        else
            t= sum(auxDiff);
        end
        disp(strcat(' * 1st time estimation inside the volume: ',num2str(durationEst)));
        disp(strcat(' * 2nd time estimation inside the volume (ts delta <= 0.3): ',num2str(t)));
        
        
        



% ====== CREATE ERRORBAR FOR Cumulative Sum of counts

p1Sum=cumsum(p1,2)';
p2Sum=cumsum(p2,2)';

p1PostCO2= p1;
p2PostCO2= p2;
filesNamePostCO2= filesName';
p1SumPostCO2= p1Sum;
p2SumPostCO2= p2Sum;

p1AIR= p1;
p2AIR= p2;
filesNameAIR= filesName';
p1SumAIR= p1Sum;
p2SumAIR= p2Sum;
%Add missing files to dataset
p1AIR= [p1AIR(1:5,:); zeros(2,length(p1AIR(1,:))); p1AIR(6:end,:)];
p2AIR= [p2AIR(1:5,:); zeros(2,length(p2AIR(1,:))); p2AIR(6:end,:)];


lg={'CO2', 'PostCo2'};
i= 1;
pp1=subplot(2,1,1);
%plot(cumsum(p1,2)');
plot(p1Sum(:,i), 'b');
pp1.YLim= [0 max(max(p1Sum(:,i)))+2000];
hold on;
plot(p1SumPostCO2(:,i), 'r');
hold off;
title('Counts in volume in position 1 over time');
xlabel('Minute #');
%xlim([1, 120]);
ylabel(' Counts');
legend(lg);
pp2=subplot(2,1,2);
plot(p2Sum(:,i), 'b');
pp2.YLim= [0 max(max(p2Sum(:,i)))+2000];
hold on;
plot(p2SumPostCO2(:,i), 'r');
hold off;
title('Counts in volume in position 2 over time');
xlabel('Minute #');
%xlim([1, 120]);
ylabel(' Counts');
legend(lg);



p1PostCO2= p1;
p2PostCO2= p2;
filesNamePostCO2= filesName';
%Check if all experiments have a value from file
p1=[p1(1:17,:); zeros(1,120); p1(18:end,:)];
p2=[p2(1:17,:); zeros(1,120); p2(18:end,:)];
p1PostCO2=[p1PostCO2(1:6,:); zeros(1,120); p1PostCO2(7:end,:); zeros(1,120)];
p2PostCO2=[p2PostCO2(1:6,:); zeros(1,120); p2PostCO2(7:end,:); zeros(1,120)];

%As experiment on 20200323 is 3 hours duration (instead of 5), withCO2 and
%postCO2 are only 1 hours long. group the time groups per couples (t and
%t+1) to make each time group == 1 minute 
p1(end,:)= horzcat(p1(end,1:2:end)+p1(end,2:2:end), zeros(1,60));
p2(end,:)= horzcat(p2(end,1:2:end)+p2(end,2:2:end), zeros(1,60));
p1PostCO2(end,:)= horzcat(p1PostCO2(end,1:2:end)+p1PostCO2(end,2:2:end), zeros(1,60));
p2PostCO2(end,:)= horzcat(p2PostCO2(end,1:2:end)+p2PostCO2(end,2:2:end), zeros(1,60));


for i=1:length(dataset)
    aa{i}=dataset(i).type;
end
uTypes= unique(aa);

i=1;
duration=60;
meanList= zeros(length(uTypes), duration);
errList= zeros(length(uTypes), duration);
meanListPostCO2= zeros(length(uTypes), duration);
errListPostCO2= zeros(length(uTypes), duration);
for ia=uTypes
    indexes= find(strcmp(aa, ia));
    %Sum the counts per position for each of the time groups
    sumPerType=p1(indexes,1:duration)+p2(indexes,1:duration);
    sumPerTypePostCO2=p1PostCO2(indexes,1:duration)+p2PostCO2(indexes,1:duration);
    if length(indexes) > 1
         %If several exp, generates the mean over their cummulative sum over time
         meanPerOdor= mean(cumsum(sumPerType,2));
         meanPerOdorPostCO2= mean(cumsum(sumPerTypePostCO2,2));
         Err= std(cumsum(sumPerType,2))/length(indexes);
         ErrPostCO2= std(cumsum(sumPerType,2))/length(indexes);
 
    else
        meanPerOdor=cumsum(sumPerType, 2);
        meanPerOdorPostCO2=cumsum(sumPerTypePostCO2, 2);
        Err= zeros(1,duration);
        ErrPostCO2= zeros(1,duration);
    end
    meanList(i,:)= meanPerOdor;
    errList(i,:)= Err;
    meanListPostCO2(i,:)= meanPerOdorPostCO2;
    errListPostCO2(i,:)= ErrPostCO2;
    i=i+1;
end


pp1=subplot(2,1,1);
errorbar(meanList(:,:)', errList')
pp1.YLim= [0 max(max(meanList))+2000];
hold on;
%errorbar(sumPerType(:,2), ErrList)
hold off;
title('Counts in volume in withCO2');
xlabel('Minute #');
%xlim([1, 120]);
ylabel(' Counts');
legend(uTypes, 'Location','eastoutside');

pp2=subplot(2,1,2);
errorbar(meanListPostCO2(:,:)', errListPostCO2')
pp2.YLim= [0 max(max(meanPerOdor, meanPerOdorPostCO2))+2000];
hold on;
%errorbar(sumPerTypePostCO2(:,2), ErrListPostCO2)
hold off;
title('Counts in volume in postCO2');
xlabel('Minute #');
%xlim([1, 120]);
ylabel(' Counts');
legend(uTypes, 'Location','eastoutside');


%Subplot each insect type per subplot
figure(gcf);
yMax= max(max(meanList));
for insTypeIndex= 1: length(uTypes)
    pp1=subplot(3,3,insTypeIndex)
    errorbar(meanList(insTypeIndex,:), errList(insTypeIndex,:))
    pp1.YLim= [0 yMax+2000];
    hold on;
    errorbar(meanListPostCO2(insTypeIndex,:), errListPostCO2(insTypeIndex,:))
    hold off;
    title(['Counts in volume over time'; strcat('Type: ',{' '},cellstr(uTypes(insTypeIndex)))]);
    xlabel('Minute #');
    ylabel(' Counts');
    legend(lg, 'Location', 'northwest');
end
save_plot_in_exp_folder(gcf, 'countsOverTime_1hr_allTypes');




%==== for 1 type of mutant only

insType='l6';
indexes= find(strcmp(aa, insType));
if length(indexes) > 1
     sumPerType=p1(indexes,:)+p2(indexes,:);
     sumPerTypePostCO2=p1PostCO2(indexes,:)+p2PostCO2(indexes,:);
     meanPerOdor= mean(cumsum(sumPerType,2));
     meanPerOdorPostCO2= mean(cumsum(sumPerTypePostCO2,2));
     ErrList= std(cumsum(sumPerType,2))/length(indexes);
     ErrListPostCO2= std(cumsum(sumPerType,2))/length(indexes);
     
else
     % If there is only 1 experiemt per this insect type, no need to generate any mean or std error
     sumPerType=p1(indexes,:)+p2(indexes,:);
     sumPerTypePostCO2=p1PostCO2(indexes,:)+p2PostCO2(indexes,:);   
     meanPerOdor= cumsum(sumPerType,2);
     meanPerOdorPostCO2= cumsum(sumPerTypePostCO2,2);
     ErrList= zeros(1,length(sumPerType(1,:)));
     ErrListPostCO2= zeros(1,length(sumPerTypePostCO2(1,:)));
end

%Plot the mean values
lg={'CO2', 'PostCO2'};
figure(gcf);
errorbar(meanPerOdor, ErrList)
%plot(cumsum(p1,2)');
ylim= [0 max(max(meanPerOdor, meanPerOdorPostCO2))+2000];
hold on;
errorbar(meanPerOdorPostCO2, ErrListPostCO2)
hold off;
title(['Counts in volume over time'; strcat('Type: ',{' '},insType)]);
xlabel('Minute #');
ylabel(' Counts');
legend(lg, 'Location', 'northwest');
%save Plot
save_plot_in_exp_folder(gcf, strcat('countsOverTime_',insType));


% ==== END for 1 type mutant




% ==== Save in a .CSV  === 

m = [234 2;671 5;735 1;264 2;346 7] 
writematrix(m,'M.csv') 
% add position value to the matrices
p1= horzcat(ones(length(p1(:,1)), 1),p1);
p2= horzcat(ones(length(p2(:,1)),1)*2,p2);
data= vertcat(p1, p2);

%If we are using experiments with 4 cues
if (nnz(p3) & nnz(p4))
    % add position value to the matrices
    p3= horzcat(ones(length(p3(:,1)), 1)*3,p1);
    p4= horzcat(ones(length(p4(:,1)),1)*4,p2);
    data= vertcat(data, p3);
    data= vertcat(data, p4);
end

path=strcat(outputPath, outputFolder);
% Columns in the csv file: Position, time grp1, time grp2, ..., time grp N
% Each row represent 1 experiment
writematrix(data,strcat(path, 'countsPerTimeGrp.csv')); 



%%%%%%%%%%%--------------------------------------------



    % -- Plot Preference indexes grouped per color --
    numOfExp= [6, 5, 3, 3, 3, 3];
    data= str2double(string(table2array(smryPIs(:,2))));
    blackRed= data(1:6);
    whiteRed= data(7:11);
    whiteGray1= [data(12), data(16), data(23)];
    whiteGray2= [data(17),data(18), data(20)];
    whiteGray3= [data(13), data(16), data(21)];
    whiteGray4= [data(14), data(19), data(22)];
    
    % Estimates they mean value for each color 
    meanPIperColor= [mean(blackRed), mean(whiteRed), mean(whiteGray1), mean(whiteGray2), mean(whiteGray3), mean(whiteGray4)];
    xAxis= { 'Black vs Red', 'White vs Red', 'White vs Gray1', 'White vs Gray2', 'White vs Gray3', 'White vs Gray4'};
    % Estimate the error for the mean values
    errList= [std(blackRed)/numOfExp(1), std(whiteRed)/numOfExp(2), std(whiteGray1)/numOfExp(3), std(whiteGray2)/numOfExp(4), std(whiteGray1)/numOfExp(5), std(whiteGray1)/numOfExp(6)];
    
    % Sort the Data & Rearrange Labels
    [sortedPIs, newIndices] = sort(meanPIperColor); % sorts in *ascending* order
    sortedLabels = xAxis(newIndices); 
    sortedErr= errList(newIndices);
    
    bar(sortedPIs);
    hold on;
    e= errorbar(sortedPIs, sortedErr, 'o');
    e.Marker = '*';
    e.LineStyle: '-';
    e.LineWidth: 50.5000;
    e.MarkerSize = 10;
    e.Color = 'red';
    e.CapSize = 15;
    hold off;
    
    set(gca,'XtickLabel', sortedLabels, 'xtick',1:length(sortedPIs));
    
    if contains(version, 'R2019')
        xtickangle(45)    
    end
    title('Preference Index towards the tested color');
    print(gcf,strcat(outputPath,  't_PIgroupedSmryPlot', '.png'),'-dpng','-r300');        % *// 300 dpi
    





data= str2double(string(table2array(smryPIs(:,2))));
bar(data);
%xAxis= string(table2array(smryPIs(:,1)));
xAxis= 'Experiments',
set(gca,'XtickLabel', xAxis);
if contains(version, 'R2019')
    xtickangle(45)    
end
title('Preference Index towards the tested color');


data= str2double(string(table2array(smryCounts(:,2:3))));
bar(data);
%xAxis= string(table2array(smryPIs(:,1)));
xAxis= 'Experiments',
set(gca,'XtickLabel', xAxis);
if contains(version, 'R2019')
    xtickangle(45)    
end
title('Counts per Position');


% Load the data in 2 separate arrays
prevCO2= str2double(string(table2array(smryTrajectories(:,2:4))));
withCO2= str2double(string(table2array(smryTrajectories(:,5:7))));
postCO2= str2double(string(table2array(smryTrajectories(:,8:10))));

%MEAN values trajectories
trjMeanList= [mean(prevCO2(1:6,1)) mean(prevCO2(7:11,1)) mean(prevCO2(12:end,1)); ...
              mean(withCO2(1:6,1)) mean(withCO2(7:11,1)) mean(withCO2(12:end,1)); ...
              mean(postCO2(1:6,1)) mean(postCO2(7:11,1)) mean(postCO2(12:end,1))];
         
%MEAN values average time of flight
avgTimeMeanList= [mean(prevCO2(1:6,2)) mean(prevCO2(7:11,2)) mean(prevCO2(12:end,2)); ...
                   mean(withCO2(1:6,2)) mean(withCO2(7:11,2)) mean(withCO2(12:end,2)); ...
                   mean(postCO2(1:6,2)) mean(postCO2(7:11,2)) mean(postCO2(12:end,2))];          
       
%MEAN values trajectories
maxTimeMeanList= [mean(prevCO2(1:6,3)) mean(prevCO2(7:11,3)) mean(prevCO2(12:end,2)); ...
                   mean(withCO2(1:6,3)) mean(withCO2(7:11,3)) mean(withCO2(12:end,2)); ...
                   mean(postCO2(1:6,3)) mean(postCO2(7:11,3)) mean(postCO2(12:end,2))];          
                      
               
%data= str2double(string(table2array(smryTrajectories(:,2:10))));

subplot(2,1,1);
bar(trjMeanList);
%xAxis= string(table2array(smryTrajectories(:,1)));
xAxis= {'PrevCO2', 'WithCO2', 'PostCO2'};
set(gca,'XtickLabel', xAxis);
if contains(version, 'R2019')
    xtickangle(45)    
end
title('mean trajectories counted (per Exp.)');
legend(expList);


subplot(2,1,2);
bar(avgTimeMeanList);
% hold on;
% bar(maxTimeMeanList);
% hold off;
xAxis= {'PrevCO2', 'WithCO2', 'PostCO2'};
set(gca,'XtickLabel', xAxis);
if contains(version, 'R2019')
    xtickangle(45)    
end
title('average trajectories duration (per Exp.)');
legend(expList);


% ---------------------
if isfile(filename)
     % File exists.
else
     % File does not exist.
end

% ---------------------------------
%code to add informantion to an output report file
% fileName= strcat('outputTest','.xlsx');
% outPath= evalin('base',outputPath); 
% outputFile= strcat(outPath, fileName);
% trajSmry= evalin('base', trajSummary);
% exp= evalin('base', outputFolder);
% writetable(trajSmry, outputFile, 'Sheet', , 'Range', 'A1');


%Load information about the output file
fileName= strcat('Mosquito_Project_Report','.xlsx');
outPath= outputPath; 
outputFile= strcat(outPath, fileName);
expSheet= outputFolder(1:length(outputFolder)-1);

%add information about # of trajectories and their duration
trajSmry= trajSummary;
header= {'date','# of trajectories', 'avg time flight' 'max time flight', '# of trajectories', 'avg time flight' 'max time flight', '# of trajectories', 'avg time flight' 'max time flight'};
auxTable= horzcat(fileNameList, num2cell(trajSmry));
output= [header; auxTable];
emptyLine= {' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};

% Add information about the # times an insect has passed inside each visual 
% cue's volume
auxTable= cell2table(cell(length(countListWithCO2(:,1)),10-(length(countListWithCO2(1,:))+2)));
if length(countListWithCO2(1,:)) == 2
    header={'date', 'Pos-1-farthest', 'Pos-2-closest', strcat('baseColorPos (',baseColor,')'), ' ', ' ', ' ', ' ', ' ', ' '};
else
   header={'date', 'Pos-1', 'Pos-2', 'Pos-3', 'Pos-4', ' ', ' ', ' ', ' ', ' '};
end
auxTable= [horzcat(fileNameList, num2cell(countListWithCO2)), num2cell(baseCueIndexList), auxTable];
output=[output; emptyLine; header; table2cell(auxTable)];

% Add information about the Preference Index and the indexes for the base
% color
auxTable= cell2table(cell(length(countListWithCO2(:,1)),10-2));
% for base workspace= evalin( 'base', 'exist(''piList'',''var'') == 1' )
if exist('piList', 'var')
    header={'date', 'PI towards color', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};
    %auxTable= zeros(length(countListWithCO2(:,1)),10-length(countListWithCO2(1,:)));
end
auxTable= [horzcat(fileNameList, num2cell(piList)), auxTable];
output=[output; emptyLine; header; table2cell(auxTable)];

%Save information in the output file
writecell(output, outputFile, 'Sheet', expSheet , 'Range', 'A1');








%find indexes in function of visual cue position
whiteCueIndexList=zeros(length(dataset(:))-1,1);
for expIndex=1:length(dataset(:))-1
    whiteCueIndex= find(strcmp(cellstr(dataset(expIndex+1).expCues([2 3],1)),'white'));
    whiteCueIndexList(expIndex)= whiteCueIndex;
end

piList= zeros(length(data(:,1)),1);
    for expIndex= 1:length(data(:,1))
        whiteCueIndex= find(strcmp(cellstr(dataset(expIndex).expCues([2 3],1)),'white'));
        totalTraj= sum(data(expIndex,:));
        pi= (data(expIndex,i)- data(expIndex,k))/totalTraj;
        piList(expIndex)= pi;
        t=i;
        i=k;
        k=t;
    end




% garbage
a= [1 4 3 4 5 6 7 8 9 4 0 9 8 7 6 5 4 3 2 4];
b= [12 3 34 5 45 6 56 7 67 8 89 0 90 1 8 8 7 43 21 4];
c= [58 3 504 4 50 3 568 56 40 3 60 20 60 2 9 1 1 11 2 4];

checkA= find(a > 2 & a < 10);
checkB= find(b > 2 & b < 10);
checkC= find(c > 2 & c < 10);

checkAB= intersect(checkA, checkB);
checkABC= intersect(checkAB, checkC);

disp(checkABC)

a= [1 4 3 4 5 6 7 8 9 4 0 9 8 7 6 5 4 3 2 4];
b= [12 3 34 5 45 6 56 7 67 8 89 0 90 1 8 8 7 43 21 4];
c= [58 3 504 4 50 3 568 56 40 3 60 20 60 2 9 1 1 11 2 4];

index =[5 6 7 8 9 10 11 12 13 14 15 16 17]
checkA2= find(a(index) > 2 & a(index) < 10);
checkB2= find(b(index) > 2 & b(index) < 10);
checkC2= find(c(index) > 2 & c(index) < 10);

[checkAB2, ia, ib]= intersect(checkA2, checkB2);
[checkABC2, ia2, ib2]= intersect(checkAB2, checkC2);

e0= [a; b; c]
e= [a(checkABC); b(checkABC); c(checkABC)]
e2= [a(checkABC2+index(1)-1); b(checkABC2+index(1)-1); c(checkABC2+index(1)-1)]






y1 = [2 2 3; 2 5 6; 2 8 9; 2 11 12];
y2 = [1 2 3; 1 2 3; 1 8 9; 1 1 1];
y3 = [3 1 3; 3 2 3; 3 3 3; 3 4 3];
y4 = [6 2 6; 6 4 6; 6 8 6; 6 10 6];
for i= 1:length(y1(:, 1))
    k= [y1(1,:); y2(1,:); y3(1,:); y4(1,:)];
    test(i)= k;
end
% z vector will ve used to separate experiments in the bar plot 
% (DISPLAY ONLY!)
z= zeros(1,length(y1(1,:)));
k= [y1(1,:);y2(1,:);y3(1,:);y4(1,:);z;y1(2,:);y2(2,:);y3(2,:);y4(2,:);z; ...
    y1(3,:);y2(3,:);y3(3,:);y4(3,:);z;y1(4,:);y2(4,:);y3(4,:);y4(4,:)];

%figure()
bar([y1; y2; y3; y4], 'stacked')

xNames=cell(1, 19);
figure();
bar(k, 'stacked')
if contains(version, 'R2019')
    xtickangle(45)    
end
for i= 1:19
    xNames{i}= num2str(i);
end
set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19]);
set(gca,'XtickLabel', xNames);



k= [p1(1,:);p2(1,:);p3(1,:);p4(1,:);p1(2,:);p2(2,:);p3(2,:);p4(2,:); ...
    p1(3,:);p2(3,:);p3(3,:);p4(3,:);p1(4,:);p2(4,:);p3(4,:);p4(4,:)];

bar(k, 'stacked')


%Plot mean values per cue position
bar([sum(y1); sum(y2); sum(y3); sum(y4)], 'stacked');

t= sum(sum(y1+y2+y3+y4));
bar([sum(y1)/t; sum(y2)/t; sum(y3)/t; sum(y4)/t], 'stacked');



l =[1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19]
mod(l,5)