
% Load instants where insect are inside a cue volume and group them by time
% after the odor was released
% Arguments:
%   - numGrps: number of time groups  to divide the data 
%	- filesPath: path to the files
%   - fileLis: List of the files to work with
% Returns:
%   - filesName: List containing the files loaded in the order they were
%                loaded
%   - pX: Matrices with the number of counts in the volume located in position X 
%         for each time group per experiment. (rows: experiment index, cols: time group)
%   - tX: Matrices with the time spent by the insects in the volume located in position X 
%         for each time group per experiment. (rows: experiment index, cols: time group)
%   - totalIDsInPX: Number of trajectories IDs that were near the cue X  

function [filesName, p1, p2, p3, p4, t1, t2, totalIDsInP1, totalIDsInP2]= load_insect_data_per_time_groups(numGrps, filesPath, filesList) 

    
    % Initialize the matrices for conts detected in each time group
    p1 = zeros(length(filesList),numGrps); 
    p2 = zeros(length(filesList),numGrps);
    p3 = zeros(length(filesList),numGrps); 
    p4 = zeros(length(filesList),numGrps);

    % Initialize the matrices for the timestamps
    t1 = zeros(length(filesList),numGrps); 
    t2 = zeros(length(filesList),numGrps);   
    
    % Initialize the matrices pour les IDs. If changing experiment type (or duration),
    % it is possible that the number of columns must be changed  
    totalIDsInP1= zeros(length(filesList),1);
    totalIDsInP2= zeros(length(filesList),1);
    
    for fileIndex= 1:length(filesList)
        fileName= filesList(fileIndex).name;
        disp(strcat(' - Working with file: ', {' '}, fileName));
        filesName(fileIndex)= {fileName(1:8)};
        % data loaded structure: [vCuePosition	objID	timeStamp NaN initialTimeStamp finalTimeStamp]
        dataFromExcel = table2array(readtable(strcat(filesPath, fileName)));
        if isempty(dataFromExcel)
            disp(strcat('   -> No counts detected in file: ',{' '}, fileName));
            % The matrcies p1..p4; t1, t2 and IDsInP1, IDsInP2 will keep its 0 values in row assigned to this empty file 
        else
            initialTS= dataFromExcel(1,4);    
            finalTS= dataFromExcel(1,5);
            grpSz= (finalTS - initialTS)/numGrps;
            grpThreshold= zeros(1, numGrps+1);
            % Load the IDs near each of the positions/cues
            indexes= find(dataFromExcel(:,1) == 1);
            tempIDs= dataFromExcel(indexes, 2)';
            tempIDs= unique(tempIDs); 
            totalIDsInP1(fileIndex,:)= length(tempIDs);
            indexes= find(dataFromExcel(:,1) == 2);
            tempIDs= dataFromExcel(indexes, 2)';
            tempIDs= unique(tempIDs); 
            totalIDsInP2(fileIndex,:)= length(tempIDs);
            
            for i=0:numGrps
                grpThreshold(i+1)= initialTS+(grpSz*i);
            end
            for grpIndex= 1:numGrps %length(grpThreshold)
                %Select the indexes from the data that verify the time group
                %conditions
                dataIndexes= find(dataFromExcel(:,3) < grpThreshold(grpIndex+1) & dataFromExcel(:,3) >= grpThreshold(grpIndex));

                p1(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 1);
                p2(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 2);
                p3(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 3);
                p4(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 4);
                % To load timestamp data (time spent by insect near thew cue
                % Sorted the subset of data by cue position
                sortedData= sortrows(dataFromExcel(dataIndexes,1:3));
                %find the first appearence of the position 2 
                split= find(sortedData(:,1) ==2,1);

                % calculate the Delta time between timestamps values 
                tsDiffP1= diff(sortedData(1:(split-1), 3));
                tsDiffP2= diff(sortedData(split:end, 3));
                % For POS-1
                % find which  sequential counts didn't happen "~=consecutively"
                out= find(tsDiffP1 > 0.05);
                t=0;
                k=1;
                if any(out)
                    for i= out
                        % the insect left the volume at a given moment
                        t= t + sum(tsDiffP1(k:(i-1)));
                        k= i+1;
                    end
                    t= t + sum(tsDiffP1(k:end));
                else
                    %the insects has been inside the volume all the time
                    t= sum(tsDiffP1);
                end
                t1(fileIndex, grpIndex)= t;
                % For POS-2
                % find which  sequential counts didn't happen "~=consecutively"
                out= find(tsDiffP2 > 0.05);
                t=0;
                k=1;
                if any(out)
                    for i= out
                        % the insect left the volume at a given moment
                        t= t + sum(tsDiffP2(k:(i-1)));
                        k= i+1;
                    end
                    t= t + sum(tsDiffP2(k:end));
                else
                    %the insects has been inside the volume all the time
                    t= sum(tsDiffP2);
                end
                t2(fileIndex, grpIndex)= t;

            end
        end
    end    
    clear tempIDs indexes
    
end