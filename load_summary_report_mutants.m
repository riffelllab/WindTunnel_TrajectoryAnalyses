% Function to load the information for all experiments stored in the report file 
% Arguments
%   - filePath: Path to file
%   - fname: Name of the file
%   - sName: Name of the sheet to use
% Returns
%   - trajSummary: table with the information regrarding the trajectories
%   - countsSummary: table with the information regarding the counts per position
%   - PISummary: table containing the Preference Index information
%   - expSummary: Table containing the type of mutant used in each experiment
function [trajSummary, countsSummary, PISummary, expSummary]= load_summary_report_mutants(filePath, fname, sName)

    %ctsHeader = {'date', 'testedColor', 'Pos-1-farthest', 'Pos-2-closest' 'baseColor-Pos'};
    %PIHeader = {'date', 'PI towards color'};
    %trajHeader= {'date', 'numOfTrajAIR', 'avgTimeFltAIR',	'maxTimeFltAIR', 'numOfTrajCO2',	'avgTimeFltCO2',	'maxTimeFltCO2',	'numOfTrajPostCO2',	'avgTimeFlt_PostCO2',	'maxTimeFltPostCO2'};
    countsSummary= table();
    PISummary= table();
    trajSummary= table();
    
    if isfile(strcat(filePath, fname))
        %Load the data from the type of experiment we want to work with
        dataFromExcel = readtable(strcat(filePath, fname), 'Sheet', sName );
        % Add the information regarding the trajectories in its own table     
        trajSummary= dataFromExcel(:,[1,4:12]);
        % Add the information regarding the counts around each cue in its own table     
        countsSummary= dataFromExcel(:,[1,13:16]);
        % Add the information regarding the preference index in its own table     
        PISummary= dataFromExcel(:,[1,17]);
        % Add the type of mutant used per experiment
        expSummary= dataFromExcel(:,1:3);
        
        %Set the correct column names for the PosCounts and PI tables
        trajSummary.Properties.VariableNames= {'date', 'numOfTrajAIR', 'avgTimeFltAIR',	'maxTimeFltAIR', 'numOfTrajCO2',	'avgTimeFltCO2',	'maxTimeFltCO2',	'numOfTrajPostCO2',	'avgTimeFlt_PostCO2',	'maxTimeFltPostCO2'};
        countsSummary.Properties.VariableNames = {'date', 'testedColor', 'Pos_1_farthest', 'Pos_2_closest', 'baseColor_Pos'};
        PISummary.Properties.VariableNames = {'date', 'PI towards color'};
        expSummary.Properties.VariableNames = {'date', 'mosqType', 'mosqGender'};

    else
        disp(strcat(' * ERROR! - This file does not exist-->', strcat(filePath, fname)));
        disp(' * FILE NOT LOADED');
    end
end
