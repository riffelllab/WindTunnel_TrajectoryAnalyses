% Function to load the information for all experiments stored in the report file 
% Arguments
%   - filePath: Path to file
%   - fName: Name of the file
%   - sName: Name of the sheet to use
% Returns
%   - trajSummary: table with the information regrarding the trajectories
%   - countsSummary: table with the information regarding the counts per position
%   - PISummary: table containing the Preference Index information
function [trajSummary, countsSummary, PISummary]= load_summary_report(filePath, fName, sName)

    %ctsHeader = {'date', 'testedColor', 'Pos-1-farthest', 'Pos-2-closest' 'baseColor-Pos'};
    %PIHeader = {'date', 'PI towards color'};
    %trajHeader= {'date', 'numOfTrajAIR', 'avgTimeFltAIR',	'maxTimeFltAIR', 'numOfTrajCO2',	'avgTimeFltCO2',	'maxTimeFltCO2',	'numOfTrajPostCO2',	'avgTimeFlt_PostCO2',	'maxTimeFltPostCO2'};
    countsSummary= table();
    PISummary= table();
    trajSummary= table();
    
    if isfile(strcat(filePath, fName))
        %fName= 'Mosquito_Project_Report.xls';

        %Load the data from the type of experiment we want to work with
        dataFromExcel = readtable(strcat(filePath, fName), 'Sheet', sName );
        % Add the information regarding the trajectories in its own table     
        trajSummary= sortrows(dataFromExcel(:,1:13));
        % Add the information regarding the counts around each cue in its own table     
        countsSummary= sortrows(dataFromExcel(:,[1,14:17]));
        % Add the information regarding the preference index in its own table     
        PISummary= sortrows(dataFromExcel(:,[1,18]));

        %Set the correct column names for the PosCounts and PI tables
        trajSummary.Properties.VariableNames= {'date', 'numOfTrajAIR', 'totalTimeFltAIR','avgTimeFltAIR', 'maxTimeFltAIR', 'numOfTrajCO2', 'totalTimeFltCO2', 'avgTimeFltCO2', 'maxTimeFltCO2', 'numOfTrajPostCO2', 'totalTimeFltPostCO2', 'avgTimeFlt_PostCO2', 'maxTimeFltPostCO2'};
        countsSummary.Properties.VariableNames = {'date', 'testedColor', 'Pos_1_farthest', 'Pos_2_closest', 'baseColor_Pos'};
        PISummary.Properties.VariableNames = {'date', 'PI_towards_color'};

    else
        disp(strcat(' * ERROR! - This file does not exist-->', strcat(filePath, fName)));
        disp(' * FILE NOT LOADED');
    end
end
