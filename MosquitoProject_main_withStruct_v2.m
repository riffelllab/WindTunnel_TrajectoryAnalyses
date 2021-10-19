
%dataset= struct with fields:
%   fileName:   Name of the file
%   type:       Type of mosquito (wt, m1, m2, m0 (mutant 1 and 2))%
%   gender:     Gender of mosquites (m, f)
%   expCues:   List of odor/visual clues and their XYZ positions
%   attr_id:    Insect ID (from flydra)
%   attr_time:  acquisition timestamp (from Flydra)   
%   attr_frame: frame number (from Flydra)
%   attr_x:     X position of the Insect (from Flydra)
%   attr_y:     Y position of the Insect (from Flydra)
%   attr_z:     Z position of the Insect (from Flydra)
%   stim:       Stimulus received at the timestamp ('AIR', CO2, postCO2)


%function datasetList= load_all_dataset()
close all;
clear variables;
format long;
%opengl('save', 'software');


% Information related to the data to analyze
%workspace = 'C:\Users\Riffelllab\Documents\Matlab\AnalysisFLYDRA\';
%workspace= 'C:\Users\dalon\Documents\MATLAB\AnalysisFLYDRA_v011920\';
workspace= 'C:\Users\dalon\Documents\GitHub\AnalysisFLYDRA_2020';

inputPath= 'D:\DataFlydra\MosquitoProject\INPUT_Data\';
outputPath= 'D:\DataFlydra\MosquitoProject\OUTPUT_Graphs\';

%subFolder= 'FLYDRA_Trajectories_sept2019\';
%subFolder= 'FLYDRA_Trajectories_oct2019\4_colors\';

%subFolder= 'FLYDRA_Trajectories_4_blacks_squared\';
%subFolder= 'FLYDRA_Trajectories_4_colors_squared\';
%subFolder= 'FLYDRA_Trajectories_4_reds_squared\';
%subFolder= 'FLYDRA_Trajectories_4_grays_squared\';

%subFolder= 'FLYDRA_Trajectories_gray_red\';
%subFolder= 'FLYDRA_Trajectories_black_red\';
subFolder= 'FLYDRA_Trajectories_white_red\';
%subFolder= 'FLYDRA_Trajectories_white_gray\';
%subFolder= 'FLYDRA_Trajectories_2_blacks\';

%subFolder= 'FLYDRA_Trajectories_mutants_lineX\';
%subFolder= 'FLYDRA_Trajectories_AnStephensi\';
%subFolder= 'FLYDRA_Trajectories_HistamineFed\';

%subFolder= 'testCalibration\';
%subFolder= 'FLYDRA_Trajectories_oct2019\black_white\';

% Load the name of the files containing the experiments datasets stored in
% subFolder
inputPath= strcat(inputPath,subFolder);
cd(inputPath)
filesList=dir('*mainbrain.h5');


% Move to the MATLAB workspace
cd(workspace);
%Create the outputFolder folder
outputFolder = load_output_folder(outputPath, subFolder);

% Variables to use in this script
loadFullDataset=true;
fps=60.0; %90.0;
% Lower time threshold (in seconds) to plot a trajectory or not
flightTimeLimit= 1.5;
%Select the Control Color value
baseColor= 'R-HUE';  %select between: black, white, R-UHE
baseColorIndexList=zeros(length(filesList),1);
testedColorList= cell(length(filesList),1);
% Test section X,Y,Z limits (start_x= -lim_x and end_x= +lim_x --> Total x distance is lim_x*2)
% These lim values are the current size of the Test Section (DO NOT CHANGE IT)
lim_x= 0.9144;
lim_y= 0.3048;
lim_z= 0.6096;

fileNameList={};
% ===========================  
% Load and clean the raw data from Flydra
for expIndex= 1:length(filesList)
        
        % Load file name to work with 
        %fileName= '20200223_110547.mainbrain.h5';
        fileName= filesList(expIndex).name;
        filePath= strcat(inputPath,fileName);
        fileNameList=vertcat(fileNameList, {fileName(1:15)});
        %disp (filePath);
        % Load all the information from the FLYDRA .h5 file    
        [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z]= load_data_from_file(filePath, loadFullDataset);
        %tempData= [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z];
        
        % If using .JSON files
        % data=jsondecode(fileread('data.txt'));
        
        
        if contains(subFolder, 'AnStephensi')
            % Load experiment settings for the an.Stephensi mosquitoes
            [cuesSetup,ts_startAIR,ts_startCO2,ts_endCO2,ts_endAIR, mType, mGender]=  load_exp_settings_anStephensi(fileName);
        else
            % Load the experiment settings (cues positions & exp timestamps)
            [cuesSetup,ts_startAIR,ts_startCO2,ts_endCO2,ts_endAIR, mType, mGender]=  load_exp_settings(fileName);            
        end
        % --- Clean the data loaded from the H5 file ---
        % Erase the X,Y and Z values that are out of the test section limits
        [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z]= erase_pos_outside_WT_v2(attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z, lim_x, lim_y, lim_z);
        % Erase the entries recorded before and after the mass flow
        % controller script is launched and the entries with timeStamp == 0.0.
        expIndexes= find(attr_time(:) >= ts_startAIR & attr_time(:) <= ts_endAIR & attr_time(:) ~=0.0);
        % Erase all insectIDs that have been detected in only 1 frame (or
        % below 0.x seconds (is it worth it to keep it?)
        
        % --- ---

        %Fill the dataset with the data from current experiment
        dataset(expIndex).fileName=    fileName;
        dataset(expIndex).type=        mType;
        dataset(expIndex).gender=       mGender;
        dataset(expIndex).expCues=    cuesSetup;
        
        dataset(expIndex).attr_id=     attr_id(expIndexes);
        dataset(expIndex).attr_time=   attr_time(expIndexes);
        dataset(expIndex).attr_frame=  attr_frame(expIndexes);
        dataset(expIndex).attr_x=      attr_x(expIndexes);
        dataset(expIndex).attr_y=      attr_y(expIndexes);
        dataset(expIndex).attr_z=      attr_z(expIndexes);

        % Estimate the timestamps values for odor stimulus ON and OFF
        % attr_time contains timestamp epochs (in seconds)
        indexPrevCO2= find(dataset(expIndex).attr_time(:) < ts_startCO2);
        indexWithCO2= find(dataset(expIndex).attr_time(:) >= ts_startCO2 & dataset(expIndex).attr_time(:) < ts_endCO2-1);
        indexPostCO2= find(dataset(expIndex).attr_time(:) >= ts_endCO2);
                
        dataset(expIndex).stim(indexPrevCO2,1)= {'AIR'};
        dataset(expIndex).stim(indexWithCO2,1)= {'CO2'};
        dataset(expIndex).stim(indexPostCO2,1)= {'postCO2'};
              
        if length(cuesSetup(:,1)) == 3
            %Load the position where the baseCue is and the color tested in experiment
            baseCueIndex=find(strcmp(cuesSetup([2,3],1), baseColor));
            baseColorIndexList(expIndex)= baseCueIndex;
            testedColorList(expIndex)= cuesSetup(find(~strcmp(cuesSetup([2,3],1), baseColor))+1,1);
        end
end

%Clear temporary variables from workspace
clear fileName mType mGender cuesSetup
clear attr_id attr_time attr_frame attr_x attr_y attr_z



% ===========================    
% Count the amount of trajectories detected while stim is prior, during and
% post CO2
trajSummary= zeros(length(dataset), 9);
for expIndex= 1:length(filesList)
    %Find the indexes fro AIR, CO2 and postCO2
    indexPrevCO2 = find(strcmp(dataset(expIndex).stim(:), 'AIR'));
    indexWithCO2 = find(strcmp(dataset(expIndex).stim(:),'CO2'));
    indexPostCO2 = find(strcmp(dataset(expIndex).stim(:),'postCO2'));
    
    % Count the total number of trajectories longer than flightTimeLimit
    % seconds. AS each row value in dataset.attr_id(:) is related to 1
    % frame we only need to check in how many frames that ID has been
    % detected
    disp(strcat(' --- FILE: ',dataset(expIndex).fileName, ' --- '));
    disp(' --- Data prev CO2 --');  
    [trajSummary(expIndex,1), trajSummary(expIndex,2), trajSummary(expIndex,3)]= count_trajectories([dataset(expIndex).attr_id(indexPrevCO2) dataset(expIndex).attr_time(indexPrevCO2)], flightTimeLimit);
    disp(' --- Data during CO2 --');
   [trajSummary(expIndex,4), trajSummary(expIndex,5), trajSummary(expIndex,6)]= count_trajectories([dataset(expIndex).attr_id(indexWithCO2), dataset(expIndex).attr_time(indexWithCO2)], flightTimeLimit);
    disp(' --- Data post CO2 --');
    [trajSummary(expIndex,7), trajSummary(expIndex,8), trajSummary(expIndex,9)]= count_trajectories([dataset(expIndex).attr_id(indexPostCO2), dataset(expIndex).attr_time(indexPostCO2)], flightTimeLimit);
    disp(' =====================================');
end



% ===========================
% Plot the heatmaps for each of the experiments. There are 3 figures for
% each experiment heatmap prior/during/post CO2being released
for expIndex= 1:length(dataset)
    %Find the indexes for AIR, CO2 and postCO2
    indexPrevCO2 = find(strcmp(dataset(expIndex).stim(:), 'AIR'));
    indexWithCO2 = find(strcmp(dataset(expIndex).stim(:),'CO2'));
    indexPostCO2 = find(strcmp(dataset(expIndex).stim(:),'postCO2'));
    
    % Plot the heatmap in for the trajectories inside the test section
    if dataset(expIndex).type(1) =='m'
        insType= 'opsin_';
        if dataset(expIndex).type(2) =='0'
            insType= strcat(insType,'1&2');    
        else
            insType= strcat(insType, dataset(expIndex).type(2));
        end
    elseif dataset(expIndex).type(1) == 'l'
        insType= strcat('line_',dataset(expIndex).type(2));
    else
        insType= 'wild_type';
    end    
    imgTitle= strcat(dataset(expIndex).fileName(1:15),'_',insType,'_',dataset(expIndex).gender);
    imgTitle= strcat(dataset(expIndex).fileName(1:15));
    plotCues= true;     %Control Flag to plot the odor/visual cues
    plot_XY_XZ_heatmaps_v5(plotCues,dataset(expIndex).attr_x(indexPrevCO2), dataset(expIndex).attr_y(indexPrevCO2), dataset(expIndex).attr_z(indexPrevCO2), dataset(expIndex).expCues, strcat(num2str(expIndex),'_1_',imgTitle,'_heatmap_prevCO2'));    
    plot_XY_XZ_heatmaps_v5(plotCues,dataset(expIndex).attr_x(indexWithCO2), dataset(expIndex).attr_y(indexWithCO2), dataset(expIndex).attr_z(indexWithCO2), dataset(expIndex).expCues, strcat(num2str(expIndex),'_2_',imgTitle,'_heatmap_withCO2'));
    plot_XY_XZ_heatmaps_v5(plotCues,dataset(expIndex).attr_x(indexPostCO2), dataset(expIndex).attr_y(indexPostCO2), dataset(expIndex).attr_z(indexPostCO2), dataset(expIndex).expCues, strcat(num2str(expIndex),'_3_',imgTitle,'_heatmap_postCO2'));   
end    





% ===========================
% Count insects (of particular type and gender) near visualClue when the CO2 is being released
% An insect ID will be counted as many times as it passes over the volume
% 
% if length(dataset(1).expCues(:,1)) == 3
%     baseCueIndexList=zeros(length(dataset(:)),1);
%     xNames= cell(1,length(dataset));
%     baseColor= 'white'; % Must be black || white || R-HUE
%     for expIndex=1:length(dataset(:))
%         baseCueIndex= find(strcmp(cellstr(dataset(expIndex).expCues([2 3],1)),baseColor));
%         baseCueIndexList(expIndex)= baseCueIndex;
%     end
% end

checkedType= 'wt';
checkedGender= 'f';
xNames= cell(1,length(dataset));
countListWithCO2= zeros(length(dataset),length(dataset(1).expCues(2:end,1)));
for expIndex= 1:length(dataset)
    % disp(strcat(' - Working with expIndex:',num2str(expIndex)));
    % Create the new file Excel file name
    disp(strcat(' * Working with group: ',{' '}, dataset(expIndex).fileName(1:15)));
    fileName= strcat(dataset(expIndex).fileName(1:8), '_countsInsideCueVol.xlsx');
    
    % Returns vector with the blobal counts per each cue position [counts in pos1, counts in pos2, ..., counts in posN]
    countListWithCO2(expIndex,:)= count_insect_in_volume_v6(strcat(outputPath, outputFolder, fileName), dataset(expIndex), flightTimeLimit, 'CO2');
    
    % Pick the experiment dates to plot later as the X axis legend
    xNames{expIndex}= strcat(dataset(expIndex).fileName(1:8),' - (',dataset(expIndex).type, ' - ', dataset(expIndex).gender, ')');
end

% Calculate (#of trajectories in visual cue) / (total # of trajectories) for all entries in countListWithCO2
plot_insect_near_vCue(countListWithCO2, subFolder, checkedType, checkedGender);
% Calculate (#of trajectories in position) / (total # of trajectories) for all entries in countListWithCO2
plot_insect_near_position(countListWithCO2, checkedType, checkedGender);
% Plots counts per position per experiment. Colors assigned to position must be added a posteriori
%plot_insect_near_position_per_exp(countListWithCO2, xNames, subFolder);


% PI FOR 2 COLORS ONLY!
% Generates Preference Index towards RED or GRAY vs a base color (white or black).
% You will need the positions where the base cue is located for each 
% dataset(X).expCues(:,1) to find indexes in function of visual cue position
% baseCueIndexList=zeros(length(dataset(:))-1,1);
% xNames= cell(1,length(dataset)-1);
% baseColor= 'black'; % Must be black or white
% for expIndex=1:length(dataset(:))-1
%     baseCueIndex= find(strcmp(cellstr(dataset(expIndex+1).expCues([2 3],1)),baseColor));
%     baseCueIndexList(expIndex)= baseCueIndex;
%     xNames{expIndex}=strcat(string(dataset(expIndex+1).expCues(2,1)),'-',string(dataset(expIndex+1).expCues(3,1)));
% end
%Identify in which position the base color is located (per each entry in dataset)
% baseCueIndexList=zeros(length(dataset(:)),1);
% xNames= cell(1,length(dataset));
% baseColor= 'white'; % Must be black or white
% testedColor= 'red';%strcmp(dataset(1).expCues(2:3,1), baseColor)
% for expIndex=1:length(dataset)
%     baseCueIndex= find(strcmp(cellstr(dataset(expIndex).expCues([2 3],1)),baseColor));
%     baseCueIndexList(expIndex)= baseCueIndex;
%     xNames{expIndex}=strcat(string(dataset(expIndex).expCues(2,1)),'-',string(dataset(expIndex).expCues(3,1)));
% end

% Generates the PreferenceIndex value per each experiment date
piList= generate_PI(countListWithCO2, baseColorIndexList);
%Change the NaN by 0
piList(isnan(piList))=0;
bar(piList)

set(gca,'XtickLabel', xNames);

if contains(version, 'R2019')
    xtickangle(45)    
end
title('Preference Index towards the tested color');



% For PI values regarding the mutants, check plot_PI_for_mutants_exp.m



% PREFERENCE INDEX for experiments grouped in function of the type of GRAY
% used in the experiment (GRAY1, GRAY2, GRAY3, GRAY4)
xNames= cell(1,4);

% GRAY 1
indexesGray1=[1 5 12];
piGrpList(1)= mean(piList(indexesGray1));
errList(1)= std(piList(indexesGray1))/length(indexesGray1);
xNames{1}='Gray1 vs White';
% GRAY 2
indexesGray2=[6 7 9];
piGrpList(2)= mean(piList(indexesGray2));
errList(2)= std(piList(indexesGray2))/length(indexesGray2);
xNames{2}='Gray2 vs White';
% GRAY 3
indexesGray3=[2 4 10];
piGrpList(3)= mean(piList(indexesGray3));
errList(3)= std(piList(indexesGray3))/length(indexesGray3);
xNames{3}='Gray3 vs White';
% GRAY 4
indexesGray4=[3 8 11];
piGrpList(4)= mean(piList(indexesGray4));
errList(4)= std(piList(indexesGray4))/length(indexesGray4);
xNames{4}='Gray4 vs White';

% Sort the Data & Rearrange Labels
[sortedPIs, newIndexes] = sort(piGrpList); % sorts in *ascending* order
sortedLabels = xNames(newIndexes); 
sortedErr= errList(newIndexes);
% Plot
bar(sortedPIs);
hold on;
e= errorbar(sortedPIs, sortedErr, 'o');
e.Marker = '*';
e.LineStyle: '-';
e.LineWidth: 50.5000;
e.MarkerSize = 10;
e.Color = 'red';
e.CapSize = 15;
hold off;

set(gca,'XtickLabel', xNames);
if contains(version, 'R2019')
    xtickangle(45)    
end
title('Preference Index towards the gray cues (grouped)');




% ===========================
% Load instants where insect are inside a cue volume and group them by time after the odor was released
%load the files list 
if (exist('outputPath', 'var') &&  exist('outputFolder', 'var'))
    % Load files
    filesPath= strcat(outputPath, outputFolder);
    cd(filesPath);
    filesList=dir('*countsInsideCueVol.xlsx');
    cd(workspace);
    % Define number of groups to generate
    numGrps=120; %=120 =60 sec per group || =4 =1800 sec per grp
    %find the counts near visual cue over time and per experiment
    [p1, p2, p3, p4]= load_ts_insect_in_volume_per_groups(numGrps, filesPath, filesList);
    % save this information in a csv file in the project folder
    save_counts_over_time_in_file(p1, p2, p3, p4)
else
    disp('  * Warning: file path and/or folder not specified in the variabes workspace ');
end
    
%BLACK_VS_RED
% as the last exp was 1hour CO2 (nstead of 2 hrs) add col i and (i+1) to
% have thesame size of seconds per group as in experiments with CO2 for 2
% hours
%pxOneHrExp= odd col vector + even col vector
p1(end,:)= horzcat(p1(end,1:2:end)+p1(end,2:2:end), zeros(1,60));
p2(end,:)= horzcat(p2(end,1:2:end)+p2(end,2:2:end), zeros(1,60));
% do the accumalite sum over the time groups (columns)
p1Sum=cumsum(p1,2)';
p2Sum=cumsum(p2,2)';

%Plot the resutls
for nameIndex= 1: length(filesList)
    lg{nameIndex}= filesList(nameIndex).name(1:8);
end
pp1=subplot(2,1,1);
%plot(cumsum(p1,2)');
plot(p1Sum);
pp1.YLim= [0 max(max(p1Sum))+2000];
hold on;
line([60, 60], pp1.YLim, 'LineWidth', 2, 'Color', 'r', 'LineStyle',':');    
hold off;
title('Counts in volume in position 1 over time');
xlabel('Minute #');
%xlim([1, 120]);
ylabel(' Counts');
legend(lg, 'Location', 'northwest');
pp2=subplot(2,1,2);
%plot(cumsum(p2,2)');
plot(p2Sum);
pp2.YLim= [0 max(max(p2Sum))+2000];
hold on;
line([60, 60], pp2.YLim, 'LineWidth', 2, 'Color', 'r', 'LineStyle',':');    
hold off;
title('Counts in volume in position 2 over time');
xlabel('Minute #');
%xlim([1, 120]);
ylabel(' Counts');
legend(lg, 'Location', 'northwest');

%Plot the  counts near visual cue over time for 1 experiment
% index=6
% plot(cumsum(p1(index,:)), 'g');
% hold on
% plot(cumsum(p2(index,:)), 'b');
% hold off
% title('Counts in volume per position over time');
% xlabel('Minute #');
% ylabel(' Counts');
% legend('Pos 1', 'Pos 2', 'Location', 'northeast');



% =============== GENERATE REPORT ===============
% generate report with all the information for this type of experiments
%countListWithCO2= zeros(length(dataset),length(dataset(1).expCues(2:end,1)));
%baseCueIndexList=zeros(length(dataset(:)),1);
%baseColor='None';
%piList=0;
generate_report(outputPath, outputFolder(1:length(outputFolder)-1), fileNameList, trajSummary, countListWithCO2, baseColorIndexList, piList)



% =============== LOAD REPORT ===============
% Load the report with all the information for all experiments with a sheet
% in the file
%[smryTrajectories, smryCounts, smryPIs, expList]= load_summary_report(outputPath, 'Mosquito_Project_Report.xls');
[smryTrajectories, smryCounts, smryPIs, expList]= load_summary_report(outputPath, 'kk.xls');




% ============ FOR WHITE-BLACK cues only ===============
% Calculate preference index (PI) for the visual cues:
% PI_n= (num. trajectories close black - num. trajectories close white)/ total traj
prefIndexPrevCO2=   zeros(length(countListPrevCO2),1);
prefIndexWithCO2= zeros(length(countListWithCO2),1);
prefIndexPostCO2= zeros(length(countListPostCO2),1);
for expIndex= 1:length(countListWithCO2)
    prefIndexPrevCO2(expIndex,1)=  ((countListPrevCO2(expIndex,2) - countListPrevCO2(expIndex,3)) / (countListPrevCO2(expIndex,2) + countListPrevCO2(expIndex,3)));
    prefIndexWithCO2(expIndex,1)=   (countListWithCO2(expIndex,2) - countListWithCO2(expIndex,3)) / (countListWithCO2(expIndex,2) + countListWithCO2(expIndex,3));
    prefIndexPostCO2(expIndex,1)=   (countListPostCO2(expIndex,2) - countListPostCO2(expIndex,3)) / (countListPostCO2(expIndex,2) + countListPostCO2(expIndex,3));
    
end

% Change the NaN values (when PI= 0/totalTraj) into 0
prefIndexPrevCO2(isnan(prefIndexPrevCO2))=0;
prefIndexWithCO2(isnan(prefIndexWithCO2))=0;
prefIndexPostCO2(isnan(prefIndexPostCO2))=0;

% Plot the diferent preference index with a different shape in funtion of
% the type of mosquitoes
separatedByType= 0; % Control flag, values: 0 or 1
figure()
hold on;
shapes=['o' 's', '+', '^', '*', 'p'];
for expIndex= 1: length(prefIndexWithCO2)
    if separatedByType
        switch dataset(expIndex).type
            case 'm0'
                scatter(1,prefIndexPrevCO2(expIndex,1), shapes(1), 'b','jitter','on', 'jitterAmount',0.15, 'DisplayName', strcat(dataset(expIndex).type, '-', dataset(expIndex).gender));
                scatter(2,prefIndexWithCO2(expIndex,1), shapes(1), 'k','jitter','on', 'jitterAmount',0.15, 'DisplayName', strcat(dataset(expIndex).type, '-', dataset(expIndex).gender));
                scatter(3,prefIndexPostCO2(expIndex,1), shapes(1), 'r','jitter','on', 'jitterAmount',0.15, 'DisplayName', strcat(dataset(expIndex).type, '-', dataset(expIndex).gender));
            case 'm1'
                scatter(1,prefIndexPrevCO2(expIndex,1), shapes(2), 'b','jitter','on', 'jitterAmount',0.15, 'DisplayName', strcat(dataset(expIndex).type, '-', dataset(expIndex).gender));
                scatter(2,prefIndexWithCO2(expIndex,1), shapes(2), 'k','jitter','on', 'jitterAmount',0.15, 'DisplayName', strcat(dataset(expIndex).type, '-', dataset(expIndex).gender));
                scatter(3,prefIndexPostCO2(expIndex,1), shapes(2), 'r','jitter','on', 'jitterAmount',0.15, 'DisplayName', strcat(dataset(expIndex).type, '-', dataset(expIndex).gender));
            case 'wt'
                scatter(1,prefIndexPrevCO2(expIndex,1), shapes(3), 'b','jitter','on', 'jitterAmount',0.15, 'DisplayName', strcat(dataset(expIndex).type, '-', dataset(expIndex).gender));
                scatter(2,prefIndexWithCO2(expIndex,1), shapes(3), 'k','jitter','on', 'jitterAmount',0.15, 'DisplayName', strcat(dataset(expIndex).type, '-', dataset(expIndex).gender));
                scatter(3,prefIndexPostCO2(expIndex,1), shapes(3), 'r','jitter','on', 'jitterAmount',0.15, 'DisplayName', strcat(dataset(expIndex).type, '-', dataset(expIndex).gender));
        end;
     else
        scatter(1,prefIndexPrevCO2(expIndex,1), shapes(expIndex), 'b','jitter','on', 'jitterAmount',0.15, 'DisplayName', strcat(dataset(expIndex).type, '-', dataset(expIndex).gender));
        scatter(2,prefIndexWithCO2(expIndex,1), shapes(expIndex), 'k','jitter','on', 'jitterAmount',0.15, 'DisplayName', strcat(dataset(expIndex).type, '-', dataset(expIndex).gender));
        scatter(3,prefIndexPostCO2(expIndex,1), shapes(expIndex), 'r','jitter','on', 'jitterAmount',0.15, 'DisplayName', strcat(dataset(expIndex).type, '-', dataset(expIndex).gender));
    end
end
hold off;
title('Preference Index prev/with/post CO2');
xlim([0 4]);
ylim([-2 2]);
%legend({'m0 - f','m0 - f','m0 - m','m1 - f','m1 - m','wt - f'});
xV= cell(1,3);
xV{1}='Prior CO2'; xV{2}='With CO2'; xV{3}='Post CO2'; 
%set(gca, 'XtickLabel', [xV{1} xV{2} xV{3}]);
set(gca, 'XTick', [1 2 3]);
set(gca, 'XtickLabel', xV);


close all;
figure()
x= linspace(1,length(countListWithCO2), length(countListWithCO2))';
yPrevCO2= (countListPrevCO2(x, 2) - countListPrevCO2(x,3))/(countListPrevCO2(x, 2) + countListPrevCO2(x,3));
yWithCO2= (countListWithCO2(x, 2) - countListWithCO2(x,3))/(countListWithCO2(x, 2) + countListWithCO2(x,3));
yPostCO2= (countListPostCO2(x, 2) - countListPostCO2(x,3))/(countListPostCO2(x, 2) + countListPostCO2(x,3));
scatter([1;1;1;1;1;1], yPrevCO2(:,1), 'b','jitter','on', 'jitterAmount',0.25);
hold on;
scatter([2;2;2;2;2;2], yWithCO2(:,6), 'g', 'jitter','on', 'jitterAmount',0.25);
scatter([3;3;3;3;3;3], yPostCO2(:,6), 'r','jitter','on', 'jitterAmount',0.25);
% plot the mean values (To use if plotting per 1 type/gender only)
%line([0.8, 1.2], [mean(yPrevCO2(:,1)), mean(yPrevCO2(:,1))], 'LineWidth', 3, 'Color', 'b', 'LineStyle',':');
%line([1.8, 2.2], [mean(yWithCO2(:,1)), mean(yWithCO2(:,1))], 'LineWidth', 3, 'Color', 'g', 'LineStyle',':');
%line([2.8, 3.2], [mean(yPostCO2(:,1)), mean(yPostCO2(:,1))], 'LineWidth', 3, 'Color', 'r', 'LineStyle',':');
hold off;
%xlim([0, 4]);
title('Preference Index prev/with/post CO2');
xV= cell(1,3);
xV{1}='Prior CO2'; xV{2}='With CO2'; xV{3}='Post CO2'; 
%set(gca, 'XtickLabel', [xV{1} xV{2} xV{3}]);
set(gca, 'XTick', [1 2 3]);
set(gca, 'XtickLabel', xV);




% ===========================
% Estimate relative activity by adding all trajectories duration and then
% dividing it by the total number of trajectories. (different estimations 
% for prior/during/post CO2) 
% WARNING: WE are dividing each of the total number of trajectories prior,
% during and post by the total number of trajectories prior CO2

activityEstList= zeros(length(dataset), 6);
% activityExpList(expIndex,:)= [trajPrevCO2, totalDurPrevCO2,trajWithCO2, totalDurWithCO2, trajPostCO2, totalDurPostCO2]
relFlightActivity= zeros(length(dataset), 3);
% relFlightActivity(expIndex,:)= [totalDurPrevCO2/totalDurPrevCO2, totalDurWithCO2/totalDurPrevCO2, totalDurPostCO2/totalDurwithCO2]
xNames= cell(1,6);
for expIndex=1:length(dataset)
    totalTraj= 0;
    totalDurationFlights= 0;
    
    %Find the indexes for AIR, CO2 and postCO2
    indexPrevCO2 = find(strcmp(dataset(expIndex).stim(:), 'AIR'));
    indexWithCO2 = find(strcmp(dataset(expIndex).stim(:),'CO2'));
    indexPostCO2 = find(strcmp(dataset(expIndex).stim(:),'postCO2'));
    
    [activityEstList(expIndex,1), activityEstList(expIndex,2)]= estimate_relative_flt_activity_v2([dataset(expIndex).attr_id(indexPrevCO2) dataset(expIndex).attr_time(indexPrevCO2)], flightTimeLimit);
    [activityEstList(expIndex,3), activityEstList(expIndex,4)]= estimate_relative_flt_activity_v2([dataset(expIndex).attr_id(indexWithCO2) dataset(expIndex).attr_time(indexWithCO2)], flightTimeLimit);
    [activityEstList(expIndex,5), activityEstList(expIndex,6)]= estimate_relative_flt_activity_v2([dataset(expIndex).attr_id(indexPostCO2) dataset(expIndex).attr_time(indexPostCO2)], flightTimeLimit);
    
    %relFlightAct(expIndex,:)= [activityEstList(expIndex,2)/activityEstList(expIndex,2), activityEstList(expIndex,4)/activityEstList(expIndex,2), activityEstList(expIndex,6)/activityEstList(expIndex,2);];
    relFlightAct(expIndex,:)= [activityEstList(expIndex,2)/activityEstList(expIndex,2), activityEstList(expIndex,4)/activityEstList(expIndex,2), activityEstList(expIndex,6)/activityEstList(expIndex,4);];

    xNames{expIndex}= strcat(dataset(expIndex).fileName(1:8),' - (',dataset(expIndex).type, ' - ', dataset(expIndex).gender, ')');
end

figure()
bp= bar(relFlightAct);
%colormap(summer(3));
set(gca,'XtickLabel', xNames);
if contains(version, 'R2019')
    xtickangle(45)    
end
l= cell(1,3);
l{1}='Prior CO2'; l{2}='With CO2'; l{3}= 'Post CO2';
legend(bp,l, 'Location', 'northeast');

%title('Estimated flight activity compared to prior-CO2 activity');
title('Estimated flight activity compared to previous activity');



%===========================
% Check the insect preference for the left or right side of the TS
checkedType= 'wt';
checkedGender= 'f';
% Change this value to 1 if you want to plot the distribution of all experiments together
allExp= 1;
megaX=[];
megaY=[];
megaZ=[];
for expIndex= 1:length(filesList(:))
    if (strcmp(dataset(expIndex).type, checkedType) && strcmp(dataset(expIndex).gender, checkedGender)) || allExp
        %Load the indexes for the XYZ positions when CO2 is released
        indexWithCO2 = find(strcmp(dataset(expIndex).stim(:),'CO2'));
        megaX= cat(1,megaX,  dataset(expIndex).attr_x(indexWithCO2));
        megaY= cat(1,megaY, dataset(expIndex).attr_y(indexWithCO2));
        megaZ= cat (1,megaZ, dataset(expIndex).attr_z(indexWithCO2));
    end;
end;
if allExp
    check_preference_in_testSection([megaX, megaY, megaZ], 'All_exp_distribution when CO2 is ON', true);
else
    check_preference_in_testSection([megaX, megaY, megaZ], strcat(checkedType,'_',checkedGender,'_exp_distribution'), true);
end;



% ===========================
% check if there are IDs detected while CO2 is released to see if there is
% a change in the path
startCO2= find(strcmp(dataset(1).stim(:),'CO2'),1);
insectID= dataset(1).attr_id(startCO2);
if length(unique(dataset(1).attr_id(startCO2-3: startCO2+3)))==1
    %the Insect ID for 4_color expriment on Oct-29 is insectID= 1203
    idIndexes= find(dataset(1).attr_id(:)== insectID);
    duration=idIndexes/fps;
    %[attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z];
    posXYZforID= [dataset(1).attr_x(idIndexes), dataset(1).attr_y(idIndexes), dataset(1).attr_y(idIndexes)];
    plot_trajectory_2D_v3(insectID, posXYZforID, 'r', dataset(1).expCues, duration, dataset(1).attr_time(idIndexes(1)));

end;






% ======================================================================
% ==================================== kk ==============================
% ======================================================================

% % ===========================
% % Load instants where insect are inside a cue volume and group them by time
% % after the odor was released
% 
% %Create input folder
% %load the outputFolder folder
% inputFolder = load_output_folder(outputPath, subFolder);
% % Load files
% filesPath= strcat(outputPath, inputFolder);
% cd(filesPath);
% filesList=dir('*countsInsideCueVol.xlsx');
% cd(workspace);
% % Define number of groups to generate
% numGrps=120; %=120 =60 sec per group || =4 =1800 sec per grp
% % Initialize the counters
% p1 = zeros(length(filesList),numGrps); 
% p2 = zeros(length(filesList),numGrps);
% p3 = zeros(length(filesList),numGrps); 
% p4 = zeros(length(filesList),numGrps);
% 
% for fileIndex= 1:length(filesList)
%     fileName= filesList(fileIndex).name;
%     disp(strcat(' - Working with file: ', {' '}, fileName));
%     
%     % data loaded structure: [vCuePosition	objID	timeStamp NaN initialTimeStamp finalTimeStamp]
%     dataFromExcel = table2array(readtable(strcat(filesPath, fileName)));
%     initialTS= dataFromExcel(1,4);    
%     finalTS= dataFromExcel(1,5);
%     %disp(strcat('  - data from excel:',num2str(dataFromExcel(10,3))));
%     %disp(strcat('  - data from excel 2:',num2str(dataFromExcel(1430,3))));
%     grpSz= (finalTS - initialTS)/numGrps;
%     grpThreshold= zeros(1, numGrps+1);
%     for i=0:numGrps
%         grpThreshold(i+1)= initialTS+(grpSz*i);
%         %grpThreshold(i+1)= (grpSz*i);
% 
%     end
%     %g=groupcounts(dataFromExcel(:,3), initialTS+grpThreshold );
%     for grpIndex= 1:numGrps %length(grpThreshold)
%         %disp(strcat('  - grpIndex:',num2str(grpIndex)));
%         %disp(strcat('  - threshold:', num2str(grpThreshold(grpIndex))));
%             
%         %Select the indexes from the data that verify the time group
%         %conditions
%         dataIndexes= find(dataFromExcel(:,3) < grpThreshold(grpIndex+1) & dataFromExcel(:,3) >= grpThreshold(grpIndex));
%         p1(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 1);
%         p2(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 2);
%         p3(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 3);
%         p4(fileIndex, grpIndex)= sum(dataFromExcel(dataIndexes,1) == 4);
%     end
% end
% 
% % plot the results of all experiments grouped
% t= sum(sum(p1+p2+p3+p4));
% bp= bar([sum(p1)/t; sum(p2)/t; sum(p3)/t; sum(p4)/t], 'stacked');
% if contains(subFolder, '4_blacks')
%     xV= cell(1,4);
%     xV{1}='Pos 1 (0.45, 0.095)'; xV{2}='Pos 2 (0.675, 0.095)'; xV{3}='Pos 3 (0.45, -0.095)';xV{4}='Pos 4 (0.675, -0.095)';
% elseif contains(subFolder, '2_blacks')
%     xV= cell(1,2);
%     xV{1}='Pos 2 (0.675, 0.095)'; xV{2}='Pos 4 (0.675, -0.095)';
% elseif contains(subFolder, 'white_red')
%     bp= bar([sum(p1)/t; sum(p2)/t;], 'stacked');
%     xV= cell(1,2);
%     xV{1}=' WHITE '; xV{2}=' RED ';
% elseif contains(subFolder, 'mutants')
%     %bp= bar([sum(p1)/t; sum(p2)/t;], 'stacked');
%     subplot(2,1,1);
%     plot(p1');
%     title('Counts for position 1 (mutants) ');
%     xlabel('Minute number');
%     ylabel('Counts');
%     subplot(2,1,2)
%     plot(p2');
%     set(gca, 'XTick', 1:(length(p1)));
%     title('Counts for position 2 (mutants) ');
%     xlabel('Minute number');
%     ylabel('Counts');
% elseif contains(subFolder, 'Histamine')
%     bp= bar([sum(p1)/t; sum(p2)/t;], 'stacked');
%     xV= cell(1,2);
%     xV{1}=' BLACK '; xV{2}=' WHITE ';
% else
%     xV= cell(1,4);
%     xV{1}='Pos 1 (0.45, 0.10)'; xV{2}='Pos 2 (0.65, 0.10)'; xV{3}='Pos 3 (0.45, -0.10)';xV{4}='Pos 4 (0.65, -0.10)';
% end
% set(gca, 'XtickLabel', xV);
% if contains(version, 'R2019')
%     xtickangle(45)    
% end
% lg= cell(1, numGrps);
% lg={'TS Grp-1', 'TS Grp-2', 'TS Grp-3', 'TS Grp-4'};
% legend(bp, lg, 'Location', 'northwest');
% title(strcat('Counts for: ',{' '}, checkedType, {' - '}, checkedGender));

% %Plot the results
% % z vector will ve used to separate experiments in the bar plot 
% % (DISPLAY ONLY!)
% z= zeros(1,length(p1(1,:)));
% 
% % Organize the data to plot it grouped by experiment day
% if contains(subFolder, '4_colors')
%     dataSorted= [p1(1,:);p2(1,:);p3(1,:);p4(1,:);z;p1(2,:);p2(2,:);p3(2,:);p4(2,:);z; ...
%                  p1(3,:);p2(3,:);p3(3,:);p4(3,:);z;p1(4,:);p2(4,:);p3(4,:);p4(4,:);z; ...
%                  p1(5,:);p2(5,:);p3(5,:);p4(5,:)];
%     
%     numberOfZs= 4;
%     xPosLabels= 1:(length(lg)*length(p1(:,1))+numberOfZs);
% elseif contains(subFolder, 'Histamine')
%     dataSorted= [p1;p2];
%     numberOfZs= 0;
% 
% else
%     dataSorted= [p1(1,:);p2(1,:);p3(1,:);p4(1,:);z;p1(2,:);p2(2,:);p3(2,:);p4(2,:);z; ...
%                  p1(3,:);p2(3,:);p3(3,:);p4(3,:);z;p1(4,:);p2(4,:);p3(4,:);p4(4,:);z;p1(5,:);p2(5,:);p3(5,:);p4(5,:)];
%     numberOfZs= 3;
%     
% 
% end
% % Plot the data
% bp= bar(dataSorted, 'stacked');
% if contains(version, 'R2019')
%     % Rotate x axis legend, if current version of Matlab allow it
%     xtickangle(45);    
% end
% % Create plot and X axis legends
% lg= cell(1, numGrps);
% xNames= cell(1, length(dataSorted));
% for i= 1:(length(lg)*length(p1(:,1))+numberOfZs) %+numberOfZs, because we are adding "numberOfZs" empty columns
%     cuePosition= mod(i, length(lg)+1);
%     if  cuePosition == 0
%         %Keep empty space where the plot will display a empty column
%         %(separation between experiments
%         xNames{i}= '';
%     else
%         % Assign values to both legends
%         lg{cuePosition}= strcat('TS Group-', num2str(cuePosition));
%         xNames{i}= strcat('POS-',num2str(cuePosition)); 
%     end
% end
% set(gca, 'XTick', 1:(length(lg)*length(p1(:,1))+numberOfZs));
% set(gca,'XtickLabel', xNames);
% grid on;
% legend(bp, lg, 'Location', 'northeast');
% title(' Counts in volume near position per experiment over time');

% % PREFERENCE INDEX for experiments grouped in function of the type of GRAY
% % used in the experiment (GRAY1, GRAY2, GRAY3, GRAY4)
% xNames= cell(1,4);
% GRAY1counts= countListWithCO2(1,1)+countListWithCO2(5,1)+countListWithCO2(12,2);
% WHITEcounts= countListWithCO2(1,2)+countListWithCO2(5,2)+countListWithCO2(12,1);
% totalTraj= GRAY1counts+WHITEcounts;
% piGrp= (GRAY1counts- WHITEcounts)/totalTraj;
% piGrpList(1)= piGrp;
% xNames{1}='Gray1 vs White';
% GRAY2counts= countListWithCO2(6,1)+countListWithCO2(7,1)+countListWithCO2(9,1);
% WHITEcounts= countListWithCO2(6,2)+countListWithCO2(7,2)+countListWithCO2(9,2);
% totalTraj= GRAY2counts+WHITEcounts;
% piGrp= (GRAY2counts- WHITEcounts)/totalTraj;
% piGrpList(2)= piGrp;
% xNames{2}='Gray2 vs White';
% GRAY3counts= countListWithCO2(2,2)+countListWithCO2(4,1)+countListWithCO2(10,2);
% WHITEcounts= countListWithCO2(2,1)+countListWithCO2(4,2)+countListWithCO2(10,1);
% totalTraj= GRAY3counts+WHITEcounts;
% piGrp= (GRAY3counts- WHITEcounts)/totalTraj;
% piGrpList(3)= piGrp;
% xNames{3}='Gray3 vs White';
% GRAY4counts= countListWithCO2(3,2)+countListWithCO2(8,2)+countListWithCO2(11,2);
% WHITEcounts= countListWithCO2(3,1)+countListWithCO2(8,1)+countListWithCO2(11,1);
% totalTraj= GRAY4counts+WHITEcounts;
% piGrp= (GRAY4counts- WHITEcounts)/totalTraj;
% piGrpList(4)= piGrp;
% xNames{4}='Gray4 vs White';
% bar(piGrpList)
% set(gca,'XtickLabel', xNames);
% if contains(version, 'R2019')
%     xtickangle(45)    
% end
% title('Preference Index towards the gray cues (grouped)');

