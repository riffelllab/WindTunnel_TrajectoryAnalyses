%Function to estimate the relative flight activity by counting the number
%of trajectories and the sum of all trajectories in a given the parts of an
%experiment
% - Arguments
%   - data: subset of data with the ID and time information for a given part (PrevCO2/WithCO2/PostCO2) of a given experiment date
%   - flighttimeLimit: minimum threshold (in seconds) to consider a set of points for an ID as an insect trajectory
% - Returns
%   - toTalTraj: the number of trajectories found inside the data subset
%   - totalDurationFlights: sum of the duration (in seconds) of all trajectories grouped in totalTraj

function [totalTraj, totalDurationFlights]= estimate_relative_flt_activity_v2(data, flightTimeLimit)
    totalTraj=0;
    totalDurationFlights=0;
    % Find all the IDs for this experiment
    uniqueID= unique(data(:,1));
    uniqueID= uniqueID';
    for objID= uniqueID()
        % Load the frames where appears the current objID
        objTime= find(data(:,1) == objID);
        %estimate the duration of the flight, by working with biggest and
        %smallest timestamps for the given objID
        duration= get_trajectory_duration(data(objTime(:),2));

        %Plot only if flight duration is bigger than 0.5 seconds
        if duration >= flightTimeLimit
            totalDurationFlights= totalDurationFlights + duration;
            totalTraj= totalTraj + 1;
        end
    end

