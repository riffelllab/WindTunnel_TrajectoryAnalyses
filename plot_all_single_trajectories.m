%data content= [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z];

function [trajDetected]= plot_all_single_trajectories(data, flightTimeLimit, fps, cluesSetup)
    %Pick the trajectory of a given insect and plot it
    uniqueID=unique(data(:,1));
    %Transpose from a column matrix to a row matrix
    uniqueID=uniqueID';
    % list containing all trajectories with duration > flightTimeLimit
    trajDetected=[];
    totalTraj=0;     %trajectories (over flightTimeLimit sec) counter
    colors= ['r', 'b','g','m','c'];
    colorCounter=1;
    for objID= uniqueID
        %Load the frames where appears the current objID
        objFrame= data(:,1) == objID;
        %estimate the duration of the flight
        framesLen=nnz(objFrame()==1);
        duration=framesLen/fps;
        
        %Plot only if flight duration is bigger than 0.5 seconds
        if duration >= flightTimeLimit
            colorCounter= colorCounter+1;
            if colorCounter== 6
                colorCounter=1;
            end;
            %Load the XYZ values for the current objID
            objXYZ=data(objFrame,[4:6]);
            % Pick the first timestamp associated to the objID
            indexTS= find(data(:,1)== objID, 1);
            % Estimate the time in seconds since the action started and
            % the objID detection
            startTime= data(indexTS,2) - data(1,2);
            % Plot trajectory in 2D (XY and XZ axis)
            plot_trajectory_2D_v3(objID, objXYZ, colors(colorCounter), cluesSetup, duration, startTime);
            totalTraj= totalTraj+1;
            trajDetected(totalTraj)= objID;
        end;
    end;
    disp(strcat(' * Total amount of trajectories over ',num2str(flightTimeLimit),' seconds: ',num2str(totalTraj)));
end