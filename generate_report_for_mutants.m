% Store all the information given as argument in a .xls file
% Arguments:
%   - outPath: Output path
%   - fileName: name of the output file
%   - expSheet: Name of the sheet assigned at this experiment type
%   - fileNameList: Column vector with the dates of each experiment
%   - typegenderList: Matrix with type of mosquito (col 1) and gender of 
%                     mosquito (col 2) per each exp (row)
%   - trajSmry: Summary containing num of trajectories, mean trajectory time 
%               duration and Max trajectory time duration (for prevC)2,
%               withCO2 and postCO2)
%   - countInVolList: Counts of insects in a volume around the visual cues
%   - baseCueIndexList: Column vector with the column index of the color
%                       base (used with data in countInVolList)
%   - testedColorList: Color used  against the base color in the experiment 
%   - piList: Column vector with the Preference Index values for each 
%             experiment date
% Output
%   - A file called Mosquito_Project_Report.xlsx with a sheet called as
%       expSheet that contains the data contained in the arguments.
%       If the files already exists, the new shee is added. If the sheets
%       already exists, it is erased witht he new information
function generate_report_for_mutants(outPath,fileName, expSheet, fileNameList, typeGenderList,trajSmry, countInVolList, baseColorIndexList, testedColorList, piList)
    
    %Build output file path
    outputFile= strcat(outPath, fileName);
    %Load sheet name and color names from Main script workspace
    %expSheet= strcat('exp_with_base_color_',evalin('base','baseColor')); %expSheet(1:length(expSheet)-1);
    baseColor= evalin('base','baseColor');
    %testedColorList= evalin('base','testedColorList');
    %outPath= evalin('base', 'outputPath'); 
    %outputFile= strcat(outPath, fileName);
    %trajSmry= evalin('base', 'trajSummary');
    %fileNameList= evalin('base', 'fileNameList');
    %countInVolList= evalin('base', 'countListWithCO2');
    %baseColoIndexList= evalin('base', 'baseColorIndexList');
    %piList= evalin('base','piList');

     %Check if the project exist
     if isfile(strcat(outputFile))
        %Load the name of the sheets existing on the file
        [t,sheetList]= xlsfinfo(outputFile);
        if any(strcmp(sheetList, expSheet))
            % if the sheet already exist in the file
            dataFromExcel = readtable(outputFile, 'sheet',expSheet);
            %Compare the exp. dates to see if the current exp. analysis are already on the Summary Report
            newEntries= setdiff(fileNameList(:), dataFromExcel.date);
            % Check if any of the current experiments has already been anlyzed in the past
            if cellfun(@isempty,newEntries)
                 disp(' * Current experiment is already in the Summary Report.');
                 disp(' * Analysis not added to the Summary Report');
            else
                % Find the indexes for newEntries experiments in fileNameList
                [newEntries,indexes1]= intersect(fileNameList(:), newEntries);               
                expData=[fileNameList(indexes1,:),  typeGenderList(indexes1,1), typeGenderList(indexes1,2), num2cell(trajSmry(indexes1,:)),testedColorList(indexes1), num2cell(countInVolList(indexes1,:)), num2cell(baseColorIndexList(indexes1)), num2cell(piList(indexes1))];
                writetable([dataFromExcel; expData], outputFile, 'sheet',expSheet);
            end
            clear dataFromExcel
        else
            %If there is no sheet named as expSheet
            expData= group_data_for_report(fileNameList, typeGenderList(:,1), typeGenderList(:,2), trajSmry, testedColorList, countInVolList, baseColorIndexList, piList, baseColor);
            writetable(expData, outputFile, 'sheet',expSheet, 'Range', 'A1');
        end  
     else
        output= group_data_for_report_mutants (fileNameList, typeGenderList, trajSmry, testedColorList, countInVolList, baseColorIndexList, piList, baseColor);
        writecell(output, outputFile, 'Sheet', expSheet , 'Range', 'A1');
     end
    end

     
     
     
%     %If we have data in the trajSmry parameter, prerare the date for storage
%     if nnz(trajSmry)     
%         header= {'date','num of trajectories', 'avg time flight' 'max time flight', 'num of trajectories', 'avg time flight' 'max time flight', 'num of trajectories', 'avg time flight' 'max time flight'};
%         auxTable= horzcat(fileNameList, num2cell(trajSmry));
%         %emptyLine= {' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};
%         output= [header; auxTable];
%     else
%         disp('  * Matrix containing trajectory information is empty')
%     end
%     
%     %If we have data ~=0 in the countInVolList parameter, prerare the date for storage
%     if nnz(countInVolList)
%         % Add information about the # times an insect has passed inside each visual 
%         % cue's volume
%         %auxTable= cell2table(cell(length(countInVolList(:,1)),10-(length(countInVolList(1,:))+2)));
%         if length(countInVolList(1,:)) == 2
%             header={'testedColor', 'Pos-1-farthest', 'Pos-2-closest', strcat('baseColorPos (',baseColor,')')};
%             %Creates a column with the value of color tested
%             tColorList= cell(length(countInVolList(:,1)),1);
%             tColorList(:)={testedColor};
%            
%         else
%            header={'date', 'Pos-1', 'Pos-2', 'Pos-3', 'Pos-4', strcat('baseColorPos (',evalin('base','baseColor'),')'), ' ', ' ', ' ', ' '};
%         end
%         auxTable= [tColorList, num2cell(countInVolList), num2cell(baseCueIndexList)];
%         auxTable= [header; auxTable];
%         output= [output, auxTable];
%         %output=[output; header; table2cell(auxTable); emptyLine];
%  %      output=[output; emptyLine; header; table2cell(auxTable)];
%     else
%         disp('  * Matrix containing counts in volume information is empty')
%     end
%    
%     %If we have data in the piList parameter, prerare the date for storage
%     if nnz(piList)
%         % Add information about the Preference Index and the indexes for the base
%         % color
%         header={'PI towards tested color'};
%         auxTable= [header; num2cell(piList)];
%         output=[output, auxTable];
%     else
%         disp('  * Matrix containing preference index information is empty')
%     end

    %Save information in the output file
%     writecell(output, outputFile, 'Sheet', expSheet , 'Range', 'A1');
% 
% end