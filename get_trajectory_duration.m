% Function to get the duration of a trajectory regarding timeStamps for the
% first and last frames assigned to an object ID

function duration = get_trajectory_duration(timeStampsObjID)
    %estimate the duration of the flight, checking operating with the
    %biggest and smallest timeStamp of the insect ID
    startTime= min(timeStampsObjID(:));
    endTime= max(timeStampsObjID(:));
    duration=endTime - startTime;
end