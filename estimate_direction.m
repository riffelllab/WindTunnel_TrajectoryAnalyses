% Estimates the direction of an insect at each timestamp. Direction is in
% degrees with respect to the positive x-axis. Returns the given cell array
% with directions dxy and dxz appended for each timestamp.
%       - Argument data is a cell array of trajectories with columns [attr_id, attr_time,
%       attr_frame, attr_x, attr_y, attr_z, stim] for an h5 file in
%       dataset.

function directions = estimate_direction(data)
    if isempty(data)
        return
    end
    flightTimeLimit= evalin('base', flightTimeLimit);
    uniqueIds= unique(cell2mat(data(:,1)));
        % direction in XY plane
        dxy= [];
        % direction in XZ plane
        dxz= [];
        % For each id, get their estimated direction at each timestamp.
        % NOTE: Direction at first timestamp for each insect is set to 0.
        for id= transpose(uniqueIds)
            idData= data(cell2mat(data(:,1)) == id, :);
            % Estimate the duration of each insect ID trajectory
            duration= get_trajectory_duration(data(idData(:),2));
            %Consider only if flight duration is bigger than flightTimeLimit seconds
            if duration >= flightTimeLimit
                x= cell2mat(idData(:,4));
                y= cell2mat(idData(:,5));
                z= cell2mat(idData(:,6));
                id_dxy= [0];
                id_dxz= [0];
                % skip first row
                for ts= 2:length(idData)
                    dx= x(ts) - x(ts-1);
                    dy= y(ts) - y(ts-1);
                    dz= z(ts) - z(ts-1);
                    angleXY= atan2d(dy,dx);
                    angleXZ= atan2d(dz,dx);
                    id_dxy= cat(1,id_dxy,angleXY);
                    id_dxz= cat(1,id_dxz,angleXZ);
                end
                dxy= cat(1,dxy,id_dxy);
                dxz= cat(1,dxz,id_dxz);
            end
        end
        %Return data with columns dxy (direction XY axis) and dxz
        %(directions XZ axis), appearing in that order
        directions= cat(2,data,num2cell(dxy),num2cell(dxz));
end
