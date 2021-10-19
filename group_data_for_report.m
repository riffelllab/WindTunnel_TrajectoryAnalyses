%Function to group all the differente matrices together as a table and
%prepare it to be saved in a file
% Arguments:
%   - fileNameList: Column vector with the dates of each experiment
%   - trajSmry: Summary containing num of trajectories, mean trajectory time 
%               duration and Max trajectory time duration (for prevC)2,
%               withCO2 and postCO2)
%   - testedColorList: Vertical list with the color used as test per experiment 
%   - countInVolList:Counts of insects in a volume around the visual cues
%   - BaseCueIndexList: Column vector with the column index of the color
%                       base (used with data in countInVolList)
%   - piList: Column vector with the Preference Index values for each 
%             experiment date
% Output
%   - output: cell array with all the data. 1 row == data from 1 exp
function output= group_data_for_report (fileNameList, trajSmry, testedColorList, countInVolList, baseCueIndexList, piList, baseColor)
    %If we have data ~=0 in the trajSmry parameter, prerare the date for storage
    if nnz(trajSmry)     
        header= {'date', 'numOfTrajectoriesAIR', 'totalTimeFlightAIR', 'avgTimeFlightAIR', 'maxTimeFlightAIR', 'numOfTrajectoriesCO2', 'totalTimeFlightCO2', 'avgTimeFlightCO2', 'maxTimeFlightCO2', 'numOfTrajectoriesPostCO2', 'totalTimeFlightPostCO2', 'avgTimeFlightPostCO2', 'maxTimeFlightPostCO2'};
        auxTable= horzcat(fileNameList, num2cell(trajSmry));
        output= [header; auxTable];
    else
        disp('  * Matrix containing trajectory information is empty')
    end
    
    %If we have data ~=0 in the countInVolList parameter, prerare the date for storage
    if nnz(countInVolList)
        % Add information about the # times an insect has passed inside each visual cue's volume
        if length(countInVolList(1,:)) == 2
            header={'testedColor', 'Pos_1_farthest', 'Pos_2_closest', strcat('baseColorPos_',baseColor)};          
        else
           header={'date', 'Pos-1', 'Pos-2', 'Pos-3', 'Pos-4', strcat('baseColorPos_',baseColor)};
        end
        auxTable= [testedColorList, num2cell(countInVolList), num2cell(baseCueIndexList)];
        auxTable= [header; auxTable];
        output= [output, auxTable];
    else
        disp('  * Matrix containing counts in volume information is empty')
    end
   
    %If we have data in the piList parameter, prerare the date for storage
    if nnz(piList)
        % Add information about the Preference Index and the indexes for the base color
        header={'PITowardsTestedColor'};
        auxTable= [header; num2cell(piList)];
        output=[output, auxTable];
    else
        disp('  * Matrix containing preference index information is empty')
    end
end

