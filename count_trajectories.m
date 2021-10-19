% Function to count the number of trajectories longer than a
% flightTimeLimit. The funciton also plot and stimation of minimum and
% maximum trajectories duration detected and an estimation of the average 
% flight duration of all trajectories detected
% Arguments:
%   - data: col vector with all the different IDs detected by Flydra
%   - fps: frames per second used by Flydra
%   - flightTimeLimit: minimun time threshold to consider a real trajectory
% Returns:
%   - totalTraj: number of different trajectories, longer that flightTime
%               limit, detected by Flydra
%   - avgTime:  Average duration (in seconds) for all trajectories detected
%               by FLydra
%   - maxTime: Longest trajectory time (in seconds)


function [totalTraj, avgTime, maxTime]= count_trajectories(data, flightTimeLimit)
    %select all the IDs from the vector
    uniqueID=unique(data(:,1));
    %Transpose from a column matrix to a row matrix
    uniqueID=uniqueID';
    % Initialice local variables
    totalTraj=0; 
    flightTimeList= nan(length(uniqueID),1);
    totalFlightTime=0;
    % For each insect ID, check if 
    for index= 1:length(uniqueID)
        objID= uniqueID(index);
        %Load the frames where appears the current objID
        objTime= find(data(:,1) == objID);
        % Estimate the duration of each insect ID trajectory
        duration= get_trajectory_duration(data(objTime(:),2));

        %Consider only if flight duration is bigger than flightTimeLimit seconds
        if duration >= flightTimeLimit
            totalTraj= totalTraj+1;
            %Keep temporary stored all trajectories duration
            flightTimeList(index)=duration;
            totalFlightTime= totalFlightTime+duration;
        end
    end
   
    % Calculate average and max trajectory duration
    avgTime= totalFlightTime/totalTraj;
    maxTime= max(flightTimeList);
    
    disp(strcat(' * Total amount of trajectories over ',num2str(flightTimeLimit),' seconds: ',num2str(totalTraj)));
    %disp(strcat(' * Average flight time for this dataset (seconds): ',num2str(totalFlightTime/length(uniqueID))));
    disp(strcat(' * Average flight time for this dataset (seconds): ',num2str(avgTime)));
    disp(strcat(' * Shortest flight duration (seconds): ',num2str(min(flightTimeList))));
    disp(strcat(' * Longest flight duration (seconds): ',num2str(max(maxTime))));
    
    disp(' ----- ');
end