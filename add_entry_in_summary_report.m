% Store all the information given as argument in a .xls file
% Arguments:
%   - outPath: Output path
%   - expSheet: Name of the sheet assigned at this experiment type
%   - fileNameList: Column vector with the dates of each experiment
%   - trajSmry: Summary containing num of trajectories, mean trajectory time 
%               duration and Max trajectory time duration (for prevC)2,
%               withCO2 and postCO2)
%   - countInVolList: Counts of insects in a volume around the visual cues
%   - baseCueIndex: Column vector with the column index of the color
%                       base (used with data in countInVolList)
%   - testedColor: Cell with the color of the tested visual cue
%   - piList: Column vector with the Preference Index values for each 
%             experiment date
% Output
%   - A file called Mosquito_Project_Report.xlsx with a sheet called as
%       expSheet that contains the data contained in the arguments.
%       If the files already exists, the new shee is added. If the sheets
%       already exists, it is erased witht he new information
function add_entry_in_summary_report(outPath, fileName, expSheet, expDate, trajSmry, countInVolList, baseCueIndex, testedColor, piList)

    %outPath= evalin('base', 'outputPath'); 
    outputFile= strcat(outPath, fileName);
    %Load sheet name and color names from Main script workspace
    %expSheet= strcat('exp_with_base_color_',evalin('base','baseColor')); %expSheet(1:length(expSheet)-1);
    baseColor= evalin('base','baseColor');
    %testedColor= evalin('base','testedColor');
    %trajSmry= evalin('base', 'trajSummary');
    %expSheet= evalin('base', 'outputFolder');
    %expSheet= expSheet(1:length(expSheet)-1);
    %expDate= {dataset.fileName(1:15)};
    %countInVolList= evalin('base', 'countListWithCO2');
    %baseCueIndex= evalin('base', 'baseCueIndex');
    %piList= evalin('base','piList');

     if isfile(outputFile)
         %Load the name of the sheets existing on the file
        [t,sheetList]= xlsfinfo(outputFile);
        if any(strcmp(sheetList, expSheet))
            % if the sheet already exist in the file
            dataFromExcel = readtable(outputFile, 'sheet',expSheet);
            % Check if the experiment has already been anlyzed in the past
            if any(find(strcmp(dataFromExcel.date, expDate)))
                 disp(' * Current experiment is already in the Summary Report.');
                 disp(' * Analysis not added to the Summary Report');
            else
                expData=[expDate, num2cell(trajSmry),testedColor, num2cell(countInVolList), num2cell(baseCueIndex), num2cell(piList)];
                %writecell([dataFromExcel; expData], outputFile, 'sheet',expSheet);
                writetable([dataFromExcel; expData], outputFile, 'sheet',expSheet);
            end
            clear dataFromExcel
        else
            %If there is no sheet named as expSheet
            expData= group_data_for_report(expDate, trajSmry, testedColor, countInVolList, baseCueIndex, piList, baseColor);
            writecell(expData, outputFile, 'sheet',expSheet, 'Range', 'A1');
        end
        clear t sheetList
      else
         expData= group_data_for_report(expDate, trajSmry, testedColor, countInVolList, baseCueIndex, piList, baseColor);
         writecell(expData, outputFile, 'sheet',expSheet, 'Range', 'A1');
     end
end

    