% =================================
% Returns trajectories that start before CO2 is active and end when CO2
% is active for a given dataEntry.
% 
% Columns of returned trajectories:
% [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z, stim]

function trajInBoth= traj_in_prevCO2_and_withCO2(dataEntry)
    % get indexes of PrevCO2 and WithCO2
    indexPrevCO2= find(strcmp(dataEntry.stim(:),'AIR'));
    indexWithCO2= find(strcmp(dataEntry.stim(:),'CO2'));
    
    % get IDs from PrevCO2 and WithCO2
    prevCO2IDs= dataEntry.attr_id(indexPrevCO2);
    withCO2IDs= dataEntry.attr_id(indexWithCO2);
    
    % get IDs that are found in both PrevCO2 and WithCO2
    idsInBoth= intersect(prevCO2IDs, withCO2IDs);
    
    % get indexes of IDs in idsInBoth
    allIdIndexes= [];
    for id= transpose(idsInBoth)
        idIndexes= find(dataEntry.attr_id(:) == id);
        % Estimate the duration of each insect ID trajectory
        duration= get_trajectory_duration(dataEntry.attr_time(idIndexes));
        %Consider only if flight duration is bigger than flightTimeLimit seconds
        if duration >= flightTimeLimit
            allIdIndexes= cat(1,allIdIndexes,idIndexes);
        end
    end
    % create trajectories for IDs in idsInBoth
    trajInBoth = [];
    trajInBoth = cat(2,trajInBoth, ...
        dataEntry.attr_id(allIdIndexes), ...
        dataEntry.attr_time(allIdIndexes), ...
        dataEntry.attr_frame(allIdIndexes), ...
        dataEntry.attr_x(allIdIndexes), ...
        dataEntry.attr_y(allIdIndexes), ...
        dataEntry.attr_z(allIdIndexes));
    trajInBoth= cat(2,num2cell(trajInBoth),dataEntry.stim(allIdIndexes));
end
    
    