numGrps= 2;
radius= 0.07;
h= 0.04;
expIndex=1;

% Create the output path to save the data
outputPath= 'D:\DataFlydra\MosquitoProject\OUTPUT_Graphs\';
if contains(subFolder, '4_grays')
    outputPath=strcat(outputPath,'exp_4_grays\');   
elseif contains(subFolder, '4_reds')
    outputPath=strcat(outputPath,'exp_4_reds\');
elseif contains(subFolder, '4_colors')
    outputPath= strcat(outputPath,'exp_4_colors\    ');
end
fileName= strcat(dataset(expIndex).fileName(1:8), '_countsInsideCueVol.xlsx');
% ===========================
% Count insects (of particular type and gender) near visualClue when the CO2 is being released
% Set allExp to 0 if you want to plot only a given type and gender.
% If you want to plot data from all experiments regardless the type and gender,
% use allExp= 1
% An insect ID will be counted as many times as it passes over the volume

function count_insec_in_volume_v5( fileName, data)
    % Table containing all the counts information. 
    % 1st row is filled with 0s and will be erased at the end of this function 
    % [positionX  OBJID, instantTimeStamp]
    insectCtrInPos= zeros(1,3);

    %Find the indexes for CO2 
    %indexPrevCO2 = find(strcmp(dataset(expIndex).stim(:), 'AIR'));
    indexWithCO2 = find(strcmp(dataset(expIndex).stim(:),'CO2'));
    % Create a subset of data when the CO2 is being released:
    % [attr_id, attr_x, attr_y, attr_z]
    objXYZwithCO2=[dataset(expIndex).attr_id(indexWithCO2), dataset(expIndex).attr_time(indexWithCO2), dataset(expIndex).attr_x(indexWithCO2), dataset(expIndex).attr_y(indexWithCO2), dataset(expIndex).attr_z(indexWithCO2)];

    %Pick the trajectory of a given insect
    uniqueID=unique(objXYZwithCO2(:,1));
    %Transpose from a column matrix to a row matrix
    uniqueID=uniqueID';

    for objID= uniqueID()
        %Load the frames where appears the current objID
    %        objID= uniqueID(index);
        objFrames= find(objXYZwithCO2(:,1) == objID);
        %calculate the duration of the flight
        duration= get_trajectory_duration(objXYZwithCO2(objFrames(:),2));
    %       flightTimeList(index)= duration;
        if duration >= flightTimeLimit
            %Load the XYZ values for the current objID
            objXYZ= objXYZwithCO2(objFrames,3:5);

            % For each visual cue used, find its [X,Y] position and the number
            % of times the insect objID has been near the cue
             for colorUsed = 2:length(dataset(expIndex).expCues(:,1))                        
                % Load XY position for 1st visual clue
                x= -lim_x  + cell2mat(dataset(expIndex).expCues(colorUsed,2));
                y= cell2mat(dataset(expIndex).expCues(colorUsed,3));
                center=[x y];
                % check if the insect is inside the given volume for each of
                % the axis separately
                checkX= find(center(1)- radius < objXYZ(:,1) & objXYZ(:,1) < center(1)+ radius);
                checkY= find(center(2)- radius < objXYZ(:,2) & objXYZ(:,2) < center(2)+ radius);
                checkZ= find(objXYZ(:,3) < h);
                % Then compare the indexes to see which indexes are inside the
                % volume in the 3 axis at the same time
                checkXY= intersect(checkX, checkY);
                checkXYZ= intersect(checkXY, checkZ);

                if ~isempty(checkXYZ) 
                    %Create a row for each of the counts inside the volume
                    % [positionX  OBJID, instantTimeStamp]
                    tableRows= [zeros(length(checkXYZ,1),1), objXYZwithCO2(checkXYZ,1), objXYZwithCO2(checkXYZ,2)];
                    tableRows(:,1)= colorUsed;

                    %Add the new rows generated to the table
                    insectCtrInPos= vertcat(insectCtrInPos, tableRows);

                end
            end
        end
    end
    % Erase the first row (it should be a roww full of zeros
    if sum(insectCtrInPos(1,:)) == 0
        insectCtrInPos(1,:)= [];
    end
    % Convert the matrix into a table
    T = array2table(insectCtrInPos, 'VariableNames',{'vCuePosition','objID','timeStamp'});
    % Wrtie data in xlsx file
    writetable(T, strcat(outputPath, fileName), 'Sheet', 1, 'Range', 'B2');
end
        
        
        
        
        