
% Load instants where insect are inside a cue volume and group them by time
% after the odor was released
% Arguments:
%   - numGrps: number of time groups  to divide the data 
%	- filesPath: path to the files
%   - fileLis: List of the files to work with
% Returns:
%   - filesName: List containing the files loaded in the order they were loaded
%   - pX: Matrices with the number of counts in the volume located in position X 
%         for each time group per experiment. (rows: experiment index, cols: time group)
function [filesName, p1, p2, p3, p4]= load_ts_insect_in_volume_per_groups(numGrps, filesPath, filesList) 

    
    % Initialize the counters
    p1 = zeros(length(filesList),numGrps); 
    p2 = zeros(length(filesList),numGrps);
    p3 = zeros(length(filesList),numGrps); 
    p4 = zeros(length(filesList),numGrps);

    for fileIndex= 1:length(filesList)
        fileName= filesList(fileIndex).name;
        disp(strcat(' - Working with file: ', {' '}, fileName));
        filesName(fileIndex)= {fileName(1:15)};
        % data loaded structure: [vCuePosition	objID	timeStamp NaN initialTimeStamp finalTimeStamp]
        dataFromExcel = table2array(readtable(strcat(filesPath, fileName)));
        if isempty(dataFromExcel)
            disp(strcat('   -> No counts detected in file: ',{' '}, fileName));
            % Add zeros to the counts inside each time group 
            p1(fileIndex,:)= 0;
            p2(fileIndex,:)= 0;
            p3(fileIndex,:)= 0;
            p4(fileIndex,:)= 0;
        else
                %Load initial and final timestamps for the odor used
                initialTS= dataFromExcel(1,4);    
                finalTS= dataFromExcel(1,5);
                if (finalTS -initialTS) > 3600 
                    % If odorChecked duration > 1hr (2hours), count visit in
                    % cues onl;y for the 1st hour
                    finalTS= initialTS + 3600;
                end
                
                %disp(strcat('  - data from excel:',num2str(dataFromExcel(10,3))));
                %disp(strcat('  - data from excel 2:',num2str(dataFromExcel(1430,3))));
                grpSz= (finalTS - initialTS)/numGrps;
                grpThreshold= zeros(1, numGrps+1);
                for i=0:numGrps
                    grpThreshold(i+1)= initialTS+(grpSz*i);
                    %grpThreshold(i+1)= (grpSz*i);

                end
                %g=groupcounts(dataFromExcel(:,3), initialTS+grpThreshold );
                for grpIndex= 1:numGrps %length(grpThreshold)
                    %disp(strcat('  - grpIndex:',num2str(grpIndex)));
                    %disp(strcat('  - threshold:', num2str(grpThreshold(grpIndex))));

                    %Select the indexes from the data that verify the time group
                    %conditions
                    dataIndexes= find(dataFromExcel(:,3) < grpThreshold(grpIndex+1) & dataFromExcel(:,3) >= grpThreshold(grpIndex));

                    % ** TODO -- change for groupCounts and test it **
                    p1(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 1);
                    p2(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 2);
                    p3(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 3);
                    p4(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 4);
                end
        end
    end
end