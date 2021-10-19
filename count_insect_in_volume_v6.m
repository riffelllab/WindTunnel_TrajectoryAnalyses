
% Count insects near visualClue when a given odor is being released. 
% An insect ID will be counted as many times as it passes over the volume. 
% Data will be stored in a .xlsx sheet and returned to main project.
% Data stored in the .xlsx file: 
%   - insectCtrInPos= [Cue position  OBJID, instantTimeStamp when objID is inside volume]
%   - [initial and final timestamps associated to the odorStim used]
%   - allCheckXYZ= [row indexes for a insect ID when its XYZ position is inside a volume ]
% Arguments:
%   - outputFile: name of the .xlsx file to create with the 
%   - dataEntry: Full dataset for a given experiment (matrix of matrices)
%   - flightTimeLimit: threshold (in seconds) to consider an insect ID 
%                      positions as a trajectory
%   - odorStim: Type of stimulus used when we want to count the insect near
%               a visual cue (AIR, CO2, POSTCO2)
% Returns:
%   - countsPerPosition: Number of insects counted in each volume
%                        sourrounding a cue in an experiment 
%                       [InsectCountedInCue1, ..., InsectCountedInCueN]
function countsPerPosition= count_insect_in_volume_v6(outputPath, outputFile, dataEntry, flightTimeLimit, odorStim)

        
%      outputFile= 'kk2.xlsx';
%      dataEntry= dataset(17);
%      flightTimeLimit= 1.5;
%      odorStim= 'CO2';

    % Define  the size of the volume
    radius= 0.07;
    h= 0.04;
    % Limit of the x axis in windTunnel
    lim_x= 0.9144;

    % Check if the analysisData folder for the exp type must be created
    if ~isfolder(strcat(outputPath))
        mkdir(outputPath);
    end
    outputFile= strcat(outputPath, outputFile);
    
    % Table containing all the counts information. 
    % 1st row is filled with 0s and will be erased at the end of this function 
    % [positionX  OBJID, instantTimeStamp]
    insectCtrInPos= zeros(1,3);
    allCheckXYZ= 0;
    %Find the indexes for the given stimuli( AIR/CO2/postCO2
    indexWithOdor = find(strcmp(dataEntry.stim(:),odorStim));

    % Find intial and final timestamps associated to this dataEntry
    initialTS= min(dataEntry.attr_time(indexWithOdor));
    finalTS= max(dataEntry.attr_time(indexWithOdor));
                
    %Pick the trajectory of a given insect
    for objID= unique(dataEntry.attr_id(indexWithOdor))' %uniqueID()
        %Load the frames where appears the current objID
        objFrames= find(dataEntry.attr_id(indexWithOdor) == objID);
        % Align the new objFrames indexes to their real position in dataEntry
        objFrames= objFrames(:)+indexWithOdor(1)-1;
        %Load the frames where appears the current objID
        %objFrames= dataEntry.attr_id(indexWithOdor) == objID; %return 0 or 1 per row regarding condition
        %objFrames= indexWithOdor(objFrames); % assign the real indexes to the values == 1
        %calculate the duration of the flight
        duration= get_trajectory_duration(dataEntry.attr_time(objFrames));

        if duration >= flightTimeLimit
            % For each visual cue used, find its [X,Y] position and the number
            % of times the insect objID has been near the cue
             for colorUsed = 2:length(dataEntry.expCues(:,1))                        
                % Load XY position for 1st visual clue
                x= -lim_x  + cell2mat(dataEntry.expCues(colorUsed,2));
                y= cell2mat(dataEntry.expCues(colorUsed,3));
                center=[x y];
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
                    indexesInDataEntry= checkXYZ(:) + objFrames(1) - 1;
                    %Create a row for each of the counts inside the volume
                    % [positionX  OBJID, instantTimeStamp]
                    tableRows= [zeros(length(indexesInDataEntry),1), dataEntry.attr_id(indexesInDataEntry), dataEntry.attr_time(indexesInDataEntry)];
                    % Add the cue order position to all the entries in
                    % tableRows (1-4 for 4 visual cues) 
                    tableRows(:,1)= (colorUsed-1);

                    %Add the new rows generated to the table
                    insectCtrInPos= vertcat(insectCtrInPos, tableRows);
                    allCheckXYZ= vertcat(allCheckXYZ, indexesInDataEntry);
                end
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
    
    
    if nnz(insectCtrInPos)
        % Convert the matrix into a table
        insectCtrTable = array2table(insectCtrInPos, 'VariableNames',{'vCuePosition','objID','timeStamp'});
        % Write data in xlsx file
        writetable(insectCtrTable, outputFile, 'Sheet', odorStim, 'Range', 'A1');
    
        % Write the initial and final timestamps the Excel File
        writetable(table(initialTS, finalTS), outputFile, 'Sheet', odorStim, 'Range', 'D1')
    
        writetable(table(allCheckXYZ), outputFile, 'Sheet', odorStim, 'Range', 'F1')
    
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
        writetable(insectCtrTable, outputFile, 'Sheet', odorStim, 'Range', 'A1');
       
        countsPerPosition= [0, 0];
    end
    
end
        
        
        
        
        