%data content= [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z];


function manage_plottings(data, flightTimeLimit, fps, groupSize, cluesSetup)
    lim_x= 0.9144;
    lim_y= 0.3048;
    lim_z= 0.6096;

    %Create the test section in the plot to add the mosquito paths
    testSectionVol= load_test_section_volumen(cluesSetup);

    %Pick the trajectory of a given insect and plot it
    uniqueID=unique(data(:,1));
    %Transpose from a column matrix to a row matrix
    uniqueID=uniqueID';
    %trajectories= [57, 53, 73, 143, 222];
    totalTraj=0;     %trajectories (over flightTimeLimit sec) counter
    colors= ['r', 'b','g','m','c'];
    colorCounter=1;
    index=0;
    if groupSize==0 | groupSize== -1
        %if groupSize value == 0, then plot the trajectories for ALL
        %inscects ID
        groupSize= length(uniqueID);
    end
    while totalTraj < groupSize
        index= index+1;
        objID= uniqueID(index);
        % Load the frames where appears the current objID
        objFrame= data(:,1) == objID;
        %estimate the duration of the flight
        framesLen=nnz(objFrame()==1);
        duration=framesLen/fps;
        
        %Plot only if flight duration is bigger than 0.5 seconds
        if duration >= flightTimeLimit
            colorCounter= colorCounter+1;
            if colorCounter== 6
                colorCounter=1;
            end
            %Load the XYZ values for the current objID
            objXYZ=data(objFrame,[4:6]);
            
            if groupSize < length(uniqueID) | groupSize== -1
                % Pick the first timestamp associated to the objID
                indexTS= find(data(:,1)== objID, 1);
                % Estimate the time in seconds since the action started and
                % the objID detection
                startTime= data(indexTS,2) - data(1,2);
                % Plot trajectory in 2D (XY and XZ axis)
                %plot_trajectory_2D(objID, objXYZ,colors(colorCounter), cluesSetup);
                plot_trajectory_2D_v3(objID, objXYZ,colors(colorCounter), cluesSetup, duration, startTime);
                % Add objID to the list of trjectories detected
                %legend('hola');
            end
            % Add the trajectory to a 3D plot containing all trajectories
            if groupSize == -1
                plot_trajectory_3D(objID, objXYZ, testSectionVol, colors(colorCounter), flightTimeLimit);
            end
            totalTraj= totalTraj +1;

        end
    end
    %disp(strcat(' * Total amount of trajectories over ',num2str(flightTimeLimit),' seconds: ',num2str(totalTraj)));
end