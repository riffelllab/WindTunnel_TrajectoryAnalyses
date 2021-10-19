%% (A.I) LOAD DATA CLEAN
% Load in DataClean only data for any rtajectory longer than flightTimeLimit
% Sort the trajectories points in function of the TS recorded (and move each trajectory at a time
%for expIndex= 1:length(dataset)
    % Select the experiment from the dataset structure to work with
    expIndex=1;
    dataClean.fileName= dataset(expIndex).fileName;
    dataClean.type= dataset(expIndex).type;
    dataClean.gender= dataset(expIndex).gender;
    dataClean.expCues= dataset(expIndex).expCues;
    dataClean.attr_id= [];
    dataClean.attr_time= [];
    dataClean.attr_frame= [];
    dataClean.attr_x= [];
    dataClean.attr_y= [];
    dataClean.attr_z= [];
    dataClean.stim= [];
    
    %load only data from trajectories during AIR and CO2
    tempIndexes= find(strcmp(dataset(expIndex).stim(:),'postCO2'));
    closeCO2= min(dataset(expIndex).attr_time(tempIndexes));
    workIndexes= find(dataset(expIndex).attr_time < closeCO2);

    dataSmry=[];
    dataset(expIndex).attr_id= dataset(expIndex).attr_id(workIndexes);
    dataset(expIndex).attr_time= dataset(expIndex).attr_time(workIndexes);
    dataset(expIndex).attr_frame= dataset(expIndex).attr_frame(workIndexes);
    dataset(expIndex).attr_x= dataset(expIndex).attr_x(workIndexes);
    dataset(expIndex).attr_y= dataset(expIndex).attr_y(workIndexes);
    dataset(expIndex).attr_z= dataset(expIndex).attr_z(workIndexes);
    dataset(expIndex).stim= dataset(expIndex).stim(workIndexes);
    
    uniqueID= unique(dataset(expIndex).attr_id);   
    for index= 1:length(uniqueID)
        objID= uniqueID(index);
        %Load the frames where appears the current objID
        objTime= find(dataset(expIndex).attr_id == objID);
        % Estimate the duration of each insect ID trajectory
        duration= get_trajectory_duration(dataset(expIndex).attr_time(objTime));
        %Consider only if flight duration is bigger than flightTimeLimit seconds
        if duration >= flightTimeLimit
            dataClean.attr_id= vertcat(dataClean.attr_id, dataset(expIndex).attr_id(objTime));
            dataClean.attr_time= vertcat(dataClean.attr_time, dataset(expIndex).attr_time(objTime));
            dataClean.attr_frame= vertcat(dataClean.attr_frame, dataset(expIndex).attr_frame(objTime));
            dataClean.attr_x= vertcat(dataClean.attr_x, dataset(expIndex).attr_x(objTime));
            dataClean.attr_y= vertcat(dataClean.attr_y, dataset(expIndex).attr_y(objTime));
            dataClean.attr_z= vertcat(dataClean.attr_z, dataset(expIndex).attr_z(objTime));
            dataClean.stim= vertcat(dataClean.stim, dataset(expIndex).stim(objTime));
        end
    end
%end
% ======================

%% (A.II) Pick 100 random trajectories in AIR and 100 random trajectories in CO2 and plot them
% Using dataClean generated above to test this approach

% Select the IDs for both samples (AIR/CO2)
clear trajAIR trajCO2
expIndex=1;
sampleSz= 20;

expIndex=1; 
tmpIndexes= find(strcmp(dataClean(expIndex).stim, 'AIR'));
listIDsAIR= unique(dataClean(expIndex).attr_id(tmpIndexes));
sampleIdxs= randsample(1:length(listIDsAIR), sampleSz);
sampleIDsAIR= listIDsAIR(sampleIdxs');

tmpIndexes= find(strcmp(dataClean(expIndex).stim, 'CO2'));
listIDsCO2= unique(dataClean(expIndex).attr_id(tmpIndexes));
sampleIdxs= randsample(1:length(listIDsCO2), sampleSz);
sampleIDsCO2= listIDsCO2(sampleIdxs');

% Group the information values for the IDs in the sample for AIR
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

% Group the XYZ values for the IDs in the sample for CO2
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

% ================================


%%

%Plot the behavior of the trajectories for the sample with AIR and with CO2
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

% =======================================
%% (A.III) Animate the trajectories in AIR.
clear trajAIR_2

maxTrajSz= max(listSizesAIR(:,2));
extraEmptyStps= (maxTrajSz - listSizesAIR(:,2));
colors=[ [1 0 0]; [0 1 0]; [0 0 1]; [0 0 0]; [1 0 1]; [1 1 0]; [0 1 1]];
for i= 1: length(listSizesAIR(:,1))   
    % add the extra empty space to its attr_x and attr_z attributes
    %select a random starting time (in the animation)
    if listSizesAIR(i,2) == maxTrajSz
        startingPnt=1;
    else
        startingPnt= randperm(extraEmptyStps(i),1)-1;
    end
    trajAIR_2(i).id= sampleIDsAIR(i);
    trajAIR_2(i).startPnt=startingPnt;
    trajAIR_2(i).plotCounter=0;
    trajAIR_2(i).attr_time= [zeros(startingPnt,1); trajAIR(i).attr_time; zeros((maxTrajSz-(startingPnt+listSizesAIR(i,2))),1)];
    trajAIR_2(i).attr_x= [zeros(startingPnt,1); trajAIR(i).attr_x; zeros((maxTrajSz-(startingPnt+listSizesAIR(i,2))),1)];
    trajAIR_2(i).attr_y= [zeros(startingPnt,1); trajAIR(i).attr_y; zeros((maxTrajSz-(startingPnt+listSizesAIR(i,2))),1)];
    trajAIR_2(i).attr_z= [zeros(startingPnt,1); trajAIR(i).attr_z; zeros((maxTrajSz-(startingPnt+listSizesAIR(i,2))),1)];
    trajAIR_2(i).attr_time2= [trajAIR(i).attr_time; zeros(extraEmptyStps(i),1)];
    trajAIR_2(i).attr_x2= [trajAIR(i).attr_x; zeros(extraEmptyStps(i),1)];
    trajAIR_2(i).attr_y2= [trajAIR(i).attr_y; zeros(extraEmptyStps(i),1)];
    trajAIR_2(i).attr_z2= [trajAIR(i).attr_z; zeros(extraEmptyStps(i),1)];
end

% "Buffer" size, number of historic lines to keep, and governs the 
% corresponding fade increments.
tailSz = 20;
%If we want plot a point each STEP points
step=1;
close all

%Initialize video
%myVideo = VideoWriter(strcat(outputPath,outputFolder,'trajectoriesWithAIR_60fps_rndSt_2.mp4')); %open video file
myVideo = VideoWriter(strcat(outputPath,outputFolder,dataClean.fileName(1:15),'_AIR_sampleSz_',mat2str(sampleSz),'_.mp4')); %open video file

myVideo.FrameRate = 60;  %can adjust this, 5 - 10 works well for me
open(myVideo)


%Start animation
figure('visible','off');

axis([-lim_x lim_x 0 lim_z])
% Array of graphics objects to store the lines. Could use a cell array.
plotObj = gobjects(1, sampleSz);
plotTails= gobjects(1, sampleSz);

% plot trajectories all starting at the same time
% Plot Odor Source
plot(-0.76, 0.20, 'ro', 'MarkerFaceColor',[0.6350 0.0780 0.1840], 'MarkerSize', 10);
hold on
% Plot Visual Cue
plot(-0.43, 0, 'o', 'MarkerFaceColor','r', 'MarkerSize', 10); 
% Plot Volume
plot([-0.500, -0.500],[0, 0.04], 'b-');
plot([-0.360, -0.360],[0, 0.04], 'b-');
plot([-0.500, -0.360],[0.04, 0.04], 'b-');
for k = 1:step:maxTrajSz
    %look for each traj ID
    for i=1:sampleSz
        if k <= listSizesAIR(i,2)
            % plot current point in animation
            currentPoint(i).x= trajAIR_2(i).attr_x2(k);
            currentPoint(i).z= trajAIR_2(i).attr_z2(k);
            delete(plotObj(i));
            plotObj(i)= plot(currentPoint(i).x, currentPoint(i).z, 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'MarkerSize', 5);
            % PLot tail
            if k > 2 && k<=tailSz+1
               currentTail(i).x= trajAIR_2(i).attr_x2(1:step:k-1);
               currentTail(i).z= trajAIR_2(i).attr_z2(1:step:k-1);
               delete(plotTails(i));
               plotTails(i)= plot(currentTail(i).x ,currentTail(i).z, '.', 'Color', [0.5 0.5 0.5]);
            elseif k>tailSz 
               currentTail(i).x= trajAIR_2(i).attr_x2(k-tailSz:step:k-1);
               currentTail(i).z= trajAIR_2(i).attr_z2(k-tailSz:step:k-1);
               delete(plotTails(i));
               plotTails(i)= plot(currentTail(i).x, currentTail(i).z, '.', 'Color', [0.5 0.5 0.5]);
              end
        else
            if (k-tailSz) <= listSizesAIR(i,2)
                delete(plotObj(i));
                currentTail(i).x= trajAIR_2(i).attr_x2(k-tailSz:step:listSizesAIR(i,2));
                currentTail(i).z= trajAIR_2(i).attr_z2(k-tailSz:step:listSizesAIR(i,2));
                delete(plotTails(i));
                plotTails(i)= plot(currentTail(i).x, currentTail(i).z, '.', 'Color', [0.5 0.5 0.5]);
            else
                delete(plotObj(i));
                delete(plotTails(i));
            end
        end
        axis([-lim_x lim_x 0 lim_z])
        % pause 1/60 second: 
        %pause(0.01667)

       % Update video
        frame = getframe(gcf); %get frame
        writeVideo(myVideo, frame);  
    end
end
delete(plotObj);
delete(plotTails);
hold off

close(myVideo);
close all



% =============================================


%% (A.IV) Animate the trajectories in CO2.
clear trajCO2_2

maxTrajSz= max(listSizesCO2(:,2));
extraEmptyStps= (maxTrajSz - listSizesCO2(:,2));
for i= 1: length(listSizesCO2(:,1))   
    % add the extra empty space to its attr_x and attr_z attributes
    %select a random starting time (in the animation)
    if listSizesCO2(i,2) == maxTrajSz
        startingPnt=1;
    else
        startingPnt= randperm(extraEmptyStps(i),1)-1;
    end
    trajCO2_2(i).id= sampleIDsCO2(i);
    trajCO2_2(i).startPnt=startingPnt;
    trajCO2_2(i).plotCounter=0;
    trajCO2_2(i).attr_time= [zeros(startingPnt,1); trajCO2(i).attr_time; zeros((maxTrajSz-(startingPnt+listSizesCO2(i,2))),1)];
    trajCO2_2(i).attr_x= [zeros(startingPnt,1); trajCO2(i).attr_x; zeros((maxTrajSz-(startingPnt+listSizesCO2(i,2))),1)];
    trajCO2_2(i).attr_y= [zeros(startingPnt,1); trajCO2(i).attr_y; zeros((maxTrajSz-(startingPnt+listSizesCO2(i,2))),1)];
    trajCO2_2(i).attr_z= [zeros(startingPnt,1); trajCO2(i).attr_z; zeros((maxTrajSz-(startingPnt+listSizesCO2(i,2))),1)];
    trajCO2_2(i).attr_time2= [trajCO2(i).attr_time; zeros(extraEmptyStps(i),1)];
    trajCO2_2(i).attr_x2= [trajCO2(i).attr_x; zeros(extraEmptyStps(i),1)];
    trajCO2_2(i).attr_y2= [trajCO2(i).attr_y; zeros(extraEmptyStps(i),1)];
    trajCO2_2(i).attr_z2= [trajCO2(i).attr_z; zeros(extraEmptyStps(i),1)];
end

% "Buffer" size, number of historic lines to keep, and governs the 
% corresponding fade increments.
tailSz = 20;
%If we want plot a point each STEP points
step=1;
close all

%Initialize video
myVideo = VideoWriter(strcat(outputPath,outputFolder,dataClean.fileName(1:15),'_CO2_sampleSz_',mat2str(sampleSz),'_.mp4')); %open video file
%myVideo = VideoWriter(strcat(outputPath,outputFolder,'kklessPoints_10_NoVisible.mp4')); %open video file

myVideo.FrameRate = 60;  %can adjust this, 5 - 10 works well for me
open(myVideo)


%Start animation
figure('visible','off');

axis([-lim_x lim_x 0 lim_z])
% Array of graphics objects to store the lines. Could use a cell array.
plotObj = gobjects(1, sampleSz);
plotTails= gobjects(1, sampleSz);

% plot trajectories all starting at the same time
% Plot Odor Source
plot(-0.76, 0.20, 'ro', 'MarkerFaceColor',[0.6350 0.0780 0.1840], 'MarkerSize', 10);
hold on
% Plot Visual Cue
plot(-0.43, 0, 'o', 'MarkerFaceColor','r', 'MarkerSize', 10); 
% Plot Volume
plot([-0.500, -0.500],[0, 0.04], 'b-');
plot([-0.360, -0.360],[0, 0.04], 'b-');
plot([-0.500, -0.360],[0.04, 0.04], 'b-');
for k = 1:step:maxTrajSz
    %look for each traj ID
    for i=1:sampleSz
        if k <= listSizesCO2(i,2)
            % plot current point in animation
            currentPoint(i).x= trajCO2_2(i).attr_x2(k);
            currentPoint(i).z= trajCO2_2(i).attr_z2(k);
            delete(plotObj(i));
            plotObj(i)= plot(currentPoint(i).x, currentPoint(i).z, 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'MarkerSize', 5);
            % PLot tail
            if k > 2 && k<=tailSz+1
               currentTail(i).x= trajCO2_2(i).attr_x2(1:step:k-1);
               currentTail(i).z= trajCO2_2(i).attr_z2(1:step:k-1);
               delete(plotTails(i));
               plotTails(i)= plot(currentTail(i).x ,currentTail(i).z, '.', 'Color', [0.5 0.5 0.5]);
            elseif k>tailSz 
               currentTail(i).x= trajCO2_2(i).attr_x2(k-tailSz:step:k-1);
               currentTail(i).z= trajCO2_2(i).attr_z2(k-tailSz:step:k-1);
               delete(plotTails(i));
               plotTails(i)= plot(currentTail(i).x, currentTail(i).z, '.', 'Color', [0.5 0.5 0.5]);
            end
        else
            if (k-tailSz) <= listSizesCO2(i,2)
                delete(plotObj(i));
                currentTail(i).x= trajCO2_2(i).attr_x2(k-tailSz:step:listSizesCO2(i,2));
                currentTail(i).z= trajCO2_2(i).attr_z2(k-tailSz:step:listSizesCO2(i,2));
                delete(plotTails(i));
                plotTails(i)= plot(currentTail(i).x, currentTail(i).z, '.', 'Color', [0.5 0.5 0.5]);
            else
                delete(plotObj(i));
                delete(plotTails(i));
            end
        end
        axis([-lim_x lim_x 0 lim_z])
        % pause 1/60 second: 
        %pause(0.01667)

       % Update video
        frame = getframe(gcf); %get frame
        writeVideo(myVideo, frame);  
    end
end
delete(plotObj);
delete(plotTails);
hold off

close(myVideo);
close all



% =============================================

%% (A.III) Animate the trajectories in AIR.
% For this, we will match the size of all trajectories to the longest one
% Each trajectory will strart at a random time
clear id_* trajAIR_2 
% Find longest trajectory.
maxTrajSz= max(listSizesAIR(:,2));
extraEmptyStps= (maxTrajSz - listSizesAIR(:,2));
colors=[ [1 0 0]; [0 1 0]; [0 0 1]; [0 0 0]; [1 0 1]; [1 1 0]; [0 1 1]];
figure();
for i= 1: length(listSizesAIR(:,1))   
    %eval(['id_' num2str(listSizesAIR(i,1)) '=animatedline("Color", colors(i,:),"LineStyle", "none","Marker", ".")']);
    eval(['id_' num2str(listSizesAIR(i,1)) '=animatedline("Color", [0.5,0.5,0.5],"LineStyle", "none","Marker", ".")']);
    % add the extra empty space to its attr_x and attr_z attributes
    %select a random starting time (in the animation)
    if listSizesAIR(i,2) == maxTrajSz
        startingPnt=1;
    else
        startingPnt= randperm(extraEmptyStps(i),1)-1;
    end
    trajAIR_2(i).id= sampleIDsAIR(i);
    trajAIR_2(i).startPnt=startingPnt;
    trajAIR_2(i).plotCounter=0;
    trajAIR_2(i).attr_time= [zeros(startingPnt,1); trajAIR(i).attr_time; zeros((maxTrajSz-(startingPnt+listSizesAIR(i,2))),1)];
    trajAIR_2(i).attr_x= [zeros(startingPnt,1); trajAIR(i).attr_x; zeros((maxTrajSz-(startingPnt+listSizesAIR(i,2))),1)];
    trajAIR_2(i).attr_y= [zeros(startingPnt,1); trajAIR(i).attr_y; zeros((maxTrajSz-(startingPnt+listSizesAIR(i,2))),1)];
    trajAIR_2(i).attr_z= [zeros(startingPnt,1); trajAIR(i).attr_z; zeros((maxTrajSz-(startingPnt+listSizesAIR(i,2))),1)];
    trajAIR_2(i).attr_time2= [trajAIR(i).attr_time; zeros(extraEmptyStps(i),1)];
    trajAIR_2(i).attr_x2= [trajAIR(i).attr_x; zeros(extraEmptyStps(i),1)];
    trajAIR_2(i).attr_y2= [trajAIR(i).attr_y; zeros(extraEmptyStps(i),1)];
    trajAIR_2(i).attr_z2= [trajAIR(i).attr_z; zeros(extraEmptyStps(i),1)];
end 

% Initialize video
%myVideo = VideoWriter(strcat(outputPath,outputFolder,'trajectoriesWithAIR_60fps_rndSt_2.mp4')); %open video file
%myVideo = VideoWriter(strcat(outputPath,outputFolder,'kk.mp4')); %open video file

%myVideo.FrameRate = 60;  %can adjust this, 5 - 10 works well for me
%open(myVideo)

%Start animation
axis([-lim_x lim_x 0 lim_z])
% plot trajectories all starting at the same time
hold on
% Plot Odor Source
plot(-0.76, 0.20, 'ro', 'MarkerFaceColor',[0.6350 0.0780 0.1840], 'MarkerSize', 10);
% Plot Visual Cue
plot(-0.43, 0, 'o', 'MarkerFaceColor','r', 'MarkerSize', 10); 
% Plot Volume
plot([-0.500, -0.500],[0, 0.04], 'b-');
plot([-0.360, -0.360],[0, 0.04], 'b-');
plot([-0.500, -0.360],[0.04, 0.04], 'b-');
for k = 1:maxTrajSz
    %look for each traj ID
    for i=1:sampleSz
        %IF we still have points to plot, add new point
        if trajAIR_2(i).plotCounter <= listSizesAIR(i,2) && k>= trajAIR_2(i).startPnt
            %load its current animatedLine and add points
            currentID= evalin('base', sprintf('id_%d',listSizesAIR(i,1))); 
            if trajAIR_2(i).plotCounter == listSizesAIR(i,2)
                currentID.Color= [0.8, 0.8, 0.8];
            else
                currentID.Color= [0.5, 0.5, 0.5];
            end
            addpoints(currentID,trajAIR_2(i).attr_x2(k),trajAIR_2(i).attr_z2(k));
            trajAIR_2(i).plotCounter= trajAIR_2(i).plotCounter +1;
%         else
%             % Change color to indicate that the traj has ended
%             currentID= evalin('base', sprintf('id_%d',listSizesAIR(i,1))); 
%             trajAIR_2(i).plotCounter= trajAIR_2(i).plotCounter +1;
%             %currentID.Color= colors(i,:).*[0.8, 0.8, 0.8];       
%             %clearpoints(currentID);
        end
 
    end

    % once we have load the currents points fo all traj IDS
    drawnow 
    pause(0.1);
    % Update video
%    frame = getframe(gcf); %get frame
%    writeVideo(myVideo, frame);  
end
hold off
disp(' * Animation ended');

%End video
%close(myVideo);



%% (A.IV) Animate the trajectories in CO2.
% For this, we will match the size of all trajectories to the longest one
% Each trajectory will strart at a random time
clear id_*

maxTrajSz= max(listSizesCO2(:,2));
extraEmptyStps= (maxTrajSz - listSizesCO2(:,2));
colors=[ [1 0 0]; [0 1 0]; [0 0 1]; [0 0 0]; [1 0 1]; [1 1 0]; [0 1 1]];

for i= 1: length(listSizesCO2(:,1))   
    eval(['id_' num2str(listSizesCO2(i,1)) '=animatedline("Color", [0.5,0.5,0.5], "LineStyle", "none","Marker", ".")']);
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
        if trajCO2_2(i).plotCounter < listSizesCO2(i,2)
            %load its current animatedLine and add points
            addpoints(currentID,trajCO2_2(i).attr_x(k),trajCO2_2(i).attr_z(k));
            trajCO2_2(i).plotCounter= trajCO2_2(i).plotCounter+1;
        elseif trajCO2_2(i).plotCounter == listSizesCO2(i,2)
        	fprintf('-- CO2 - Ending trajID: %d with length %d (%d)\n',listSizesCO2(i,1), trajCO2_2(i).plotCounter, listSizesCO2(i,2));
        else
            % Change color to indicate that the traj has ended
            currentID.Color= [0.8, 0.8, 0.8];       
            %currentID.Color= colors(i,:).*[0.8, 0.8, 0.8];    
        end
        if mod(k,100)==0
            fprintf('CO2 - iteration %d for trajID: %d with Color index %d \n',k, listSizesCO2(i,1), i);
        end
    end
    % once we have load the currents points fo all traj IDS
    drawnow
    
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);  
end
hold off
    
close(myVideo);






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
        %zThreshold= lim_z/2;
        % duration >= flightTimeLimit && all(dataset(expIndex).attr_z(workIndexes(objTime)) <= zThreshold)
        
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
% ======================


%%
%%ANIMATION TEST


clear id_*;

szKk=100;
x=1:szKk;
y=szKk:-1:1;
z=ones(1,szKk);
o= x*2;
test(1).data=x;
test(2).data=y;
test(3).data=z;
test(4).data=o;
color=['r','g','b','k','m','y','c'];
color2=[ [1 0 0]; [0 1 0]; [0 0 1]; [0 0 0]; [1 0 1]; [1 1 0]; [0 1 1]];
 for i=1:4
     kk(i).id=i;
     kk(i).x= test(1).data;
     kk(i).y= test(i).data;
     
     eval(['id_' num2str(i) '=animatedline("Color", color2(i,:), "LineStyle", "none","Marker", ".")']);
 end
  
 for k = 1:szKk
    for i=1:4
        %look fo each traj ID
        currentID= evalin('base', sprintf('id_%d',i)); 
        if k < szKk
            %load its current animatedLine and add points
            addpoints(currentID,kk(i).x(k),kk(i).y(k));
        else
            % Change color to indicate that the traj has ended
            currentID.Color= color2(i,:).*[0.8, 0.8, 0.8];         end
    end
    % once we have load the currents points fo all traj IDS
    drawnow 
    
    pause(0.01)

 end


%%
%data of the problem 
t = 0:.01:2*pi;
x = cos(2*t).*(cos(t).^2);
y = sin(2*t).*(sin(t).^2);
comet(x,y, 0.0005);
 
 
for i = 1:length(x)
    plot(x(i),y(i),'Or')
    axis([min(x) max(x) min(y) max(y)]) ;
    drawnow
    pause(0.01) ; % if you want to slow the plot
end


t = 0:pi/50:4*pi;
x = -sin(t) - sin(t/2);
y = -cos(t) + cos(t/2);
p = 0.05;
comet(x,y,p)

axis([-lim_x lim_x 0 lim_z])
figure()
comet(trajAIR_2(1).attr_x, trajAIR(1).attr_z)
hold on
comet(trajAIR_2(2).attr_x, trajAIR(2).attr_z)
comet(trajAIR_2(3).attr_x, trajAIR(3).attr_z)
hold off

%%

figure;
axes('position',[0 0 1 1]);
plot1 = scatter(lon(1),lat(1),30,concentration(1),'.');
xlim([lonMin lonMax]);
ylim([latMin latMax]);
set(gca,'Color','none');
set(gca,'CLim',[0, 1E-4]);

for k = 2:length(lat) 
     plot1.XData = lon(k); 
     plot1.YData = lat(k); 
     plot1.CData = concentration(k); 
     % pause 2/10 second: 
     pause(0.2)
end
  
for t= 1:20
    for tt=max(1,t-nFade):t
        tVal=  max(0, 1 - (t-tt)/nFade);
        fprintf('t iter: %d -- tt iter: %d, tVal: %d',t, tt, tVal);
        disp('...')
    end
end





%%



% "Buffer" size, number of historic lines to keep, and governs the 
% corresponding fade increments.
tailSz = 5;

%Start animation
figure();
axis([-lim_x lim_x 0 lim_z])
% Array of graphics objects to store the lines. Could use a cell array.
plotObj = gobjects(1, sampleSz);
plotTails= gobjects(1, sampleSz);

for k = 1:maxTrajSz
    %look for each traj ID
    for i=1%:sampleSz
        %plotObj(k)=plot([trajAIR_2(i).attr_x2(1), trajAIR_2(i).attr_x2(k)],[trajAIR_2(i).attr_z2(1), trajAIR_2(i).attr_z2(k)],'ok');
        %plotObj(i)= scatter([trajAIR_2(i).attr_x2(1),trajAIR_2(i).attr_x2(k)], [trajAIR_2(i).attr_z2(1), trajAIR_2(i).attr_z2(k)],'.', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
        currentPoint(i).x= trajAIR_2(i).attr_x2(k);
        currentPoint(i).z= trajAIR_2(i).attr_z2(k);
        delete(plotObj(i));
        plotObj(i)= plot(currentPoint(i).x, currentPoint(i).z, 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'MarkerSize', 5);
        if k > 2 && k<=tailSz+1
           %plotObj(i)= scatter([trajAIR_2(i).attr_x2(1),trajAIR_2(i).attr_x2(k-1)], [trajAIR_2(i).attr_z2(1), trajAIR_2(i).attr_z2(k-1)],'.', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
           %[plotObj(i).MarkerFaceAlpha]= 0.5;
           currentTail(i).x= trajAIR_2(1).attr_x2(1:k-1);
           currentTail(i).z= trajAIR_2(1).attr_z2(1:k-1);
           delete(plotTails(i));
           plotTails(i)= plot(currentTail(i).x ,currentTail(i).z, '.', 'Color', [0.5 0.5 0.5]);
           %plot(trajAIR_2(i).attr_x2(k), trajAIR_2(i).attr_z2(k), '.', 'Color', [0 0 0]);
        elseif k>tailSz && k >= listSizesAIR(i,2)
           currentTail(i).x= trajAIR_2(i).attr_x2(k-tailSz:k-1);
           currentTail(i).z= trajAIR_2(i).attr_z2(k-tailSz:k-1);
           delete(plotTails(i));
           plotTails(i)= plot(currentTail(i).x, currentTail(i).z, '.', 'Color', [0.5 0.5 0.5]);
           %plotObj(i)= scatter([trajAIR_2(i).attr_x2(1),trajAIR_2(i).attr_x2(k-tailSz)], [trajAIR_2(i).attr_z2(1), trajAIR_2(i).attr_z2(k-tailSz)],'.', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
           %[plotObj(i).MarkerFaceAlpha]= 0;
           %plotHandler= plot([trajAIR_2(i).attr_x2(1), trajAIR_2(i).attr_x2(k-tailSz)],[trajAIR_2(i).attr_z2(1), trajAIR_2(i).attr_z2(k-tailSz)]);
           %delete(hc);
           %plot([trajAIR_2(i).attr_x2(k-tailSz+1), trajAIR_2(i).attr_x2(k-1)],[trajAIR_2(i).attr_z2(k-tailSz+1), trajAIR_2(i).attr_z2(k-1)], '.', 'Color', [0.5 0.5 0.5]);
           %plot(trajAIR_2(i).attr_x2(k), trajAIR_2(i).attr_z2(k), '.', 'Color', [0 0 0]);
        end
        axis([-lim_x lim_x 0 lim_z])
             % pause 2/10 second: 
        pause(0.01)
    end
end



