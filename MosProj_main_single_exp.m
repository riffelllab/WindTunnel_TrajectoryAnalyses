
%Script to load 1 single .h5 file and analyze it 
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

%% ========== USER PARAMETERS ==========
close all;
clear variables;
format short;
%opengl('save', 'software');
%fileID= '20200303_110206';      %Gray9.5
%fileID= '20200429_081211';      %Black
fileID= '20200911_121039';       %mutants
fileName= strcat(fileID,'.mainbrain.h5');
% Define the base color to use in the analysis
baseColor= 'white'; % Must be black || white || R-HUE
% Flag to plot the visual cues in the heatmaps and other plots or not
plotDataFlag= true;

% Information related to the data to analyze
expDep= 'expDeprecated\';
workspace= 'C:\Users\dalon\Documents\GitHub\AnalysisFLYDRA_2020';
inputPath= 'D:\DataFlydra\MosquitoProject\INPUT_Data\';
outputPath= 'D:\DataFlydra\MosquitoProject\OUTPUT_Graphs\';


%subFolder= 'FLYDRA_Trajectories_gray_red\';
%subFolder= 'FLYDRA_Trajectories_black_red\';
%subFolder= 'FLYDRA_Trajectories_white_red\';
%subFolder= 'FLYDRA_Trajectories_white_black\';
%subFolder= 'FLYDRA_Trajectories_white_blue\';
%subFolder= 'FLYDRA_Trajectories_2_blacks\';
subFolder= 'FLYDRA_Trajectories_mutants_lineX\';
%subFolder= 'FLYDRA_Trajectories_AnStephensi\';
%subFolder= 'FLYDRA_Trajectories_CxQuinquefasciatus\';
%subFolder= 'FLYDRA_Trajectories_HistamineFed\';
%subFolder= 'FLYDRA_Trajectories_white_green\';
%subFolder= 'FLYDRA_Trajectories_white_purple\';
%subFolder= 'FLYDRA_Trajectories_white_gray_4_0\';
%subFolder= 'FLYDRA_Trajectories_white_gray\';

%subFolder= 'testCalibration\';
%subFolder= 'FLYDRA_Trajectories_oct2019\black_white\';

% Variables to use in this script
loadFullDataset=true;
fps=60.0; %90.0;
% Lower time threshold (in seconds) to plot a trajectory or not
flightTimeLimit= 1.5;

% Test section X,Y,Z limits (start_x= -lim_x and end_x= +lim_x --> Total x distance is lim_x*2)
% These lim values are the current size of the Test Section (DO NOT CHANGE IT)
lim_x= 0.9144;
lim_y= 0.3048;
lim_z= 0.6096;

%Create the outputFolder folder
outputFolder = load_output_folder(outputPath, subFolder);
% ===============================

%% ------------------------------------------------------------------------
% METHODS TO GENERATE SUMMARIES/FILES FROM RAW DATA
% -------------------------------------------------------------------------

% Specify the  path of the file to work with 
filePath= strcat(inputPath, subFolder, fileName);
% Load all the information from the FLYDRA .h5 file    
[attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z]= load_data_from_file(filePath, loadFullDataset);
%tempData= [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z];

if contains(subFolder, 'AnStephensi')
    % Load experiment settings for the an.Stephensi mosquitoes
    [cuesSetup,ts_startAIR,ts_startCO2,ts_endCO2,ts_endAIR, mType, mGender]=  load_exp_settings_anStephensi(fileName);
elseif contains(subFolder, 'CxQuinquefasciatus')
    [cuesSetup,ts_startAIR,ts_startCO2,ts_endCO2,ts_endAIR, mType, mGender]=  load_exp_settings_cxQuinquefasciatus(fileName);
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
dataset.fileName=    fileName;
dataset.type=        mType;
dataset.gender=       mGender;
dataset.expCues=    cuesSetup;

dataset.attr_id=     attr_id(expIndexes);
dataset.attr_time=   attr_time(expIndexes);
dataset.attr_frame=  attr_frame(expIndexes);
dataset.attr_x=      attr_x(expIndexes);
dataset.attr_y=      attr_y(expIndexes);
dataset.attr_z=      attr_z(expIndexes);

% Estimate the timestamps values for odor stimulus ON and OFF
% attr_time contains timestamp epochs (in seconds)
indexPrevCO2= find(dataset.attr_time(:) < ts_startCO2);
indexWithCO2= find(dataset.attr_time(:) >= ts_startCO2 & dataset.attr_time(:) < ts_endCO2);
indexPostCO2= find(dataset.attr_time(:) >= ts_endCO2);

% Define which odor stime was active at each timestamp in teh dataset
dataset.stim(indexPrevCO2,1)= {'AIR'};
dataset.stim(indexWithCO2,1)= {'CO2'};
dataset.stim(indexPostCO2,1)= {'postCO2'};

%Load the position where the baseCue is and the color tested in experiment
baseColorIndex=find(strcmp(cuesSetup([2,3],1), baseColor));
testedColor= cuesSetup(find(~strcmp(cuesSetup([2,3],1), baseColor))+1,1);

%Clear temporary variables from workspace
clear fileName loadFullDataset indexPrevCO2 indexWithCO2 indexPostCO2
clear attr_id attr_time attr_frame attr_x attr_y attr_z
clear  cuesSetup ts_startAIR ts_startCO2 ts_endCO2 ts_endAIR mType mGender


%% ===========================
% Estimate relative activity by adding all trajectories duration and then
% dividing it by the total number of trajectories. (different estimations 
% for prior/during/post CO2) 
% WARNING: WE are dividing each of the total number of trajectories prior,
% during and post by the total number of trajectories prior CO2
activityEstList= zeros(length(dataset), 6);
relFlightActDuration= zeros(length(dataset), 3);
xNames= cell(1,length(dataset));
for expIndex=1:length(dataset)
    totalTraj= 0;
    totalDurationFlights= 0;
    
    %Find the indexes for AIR, CO2 and postCO2
    indexPrevCO2 = find(strcmp(dataset(expIndex).stim(:), 'AIR'));
    indexWithCO2 = find(strcmp(dataset(expIndex).stim(:),'CO2'));
    indexPostCO2 = find(strcmp(dataset(expIndex).stim(:),'postCO2'));
    % columns 1/3/5 has the # of trajectories detected per PrevCO2/WithCO2/PostCO2
    % and 2/4/6 contains the sum of these trajectories duration (in seconds)
    [activityEstList(expIndex,1), activityEstList(expIndex,2)]= estimate_relative_flt_activity_v2([dataset(expIndex).attr_id(indexPrevCO2) dataset(expIndex).attr_time(indexPrevCO2)], flightTimeLimit);
    [activityEstList(expIndex,3), activityEstList(expIndex,4)]= estimate_relative_flt_activity_v2([dataset(expIndex).attr_id(indexWithCO2) dataset(expIndex).attr_time(indexWithCO2)], flightTimeLimit);
    [activityEstList(expIndex,5), activityEstList(expIndex,6)]= estimate_relative_flt_activity_v2([dataset(expIndex).attr_id(indexPostCO2) dataset(expIndex).attr_time(indexPostCO2)], flightTimeLimit);
    
    %relFlightActDuration(expIndex,:)= [activityEstList(expIndex,2)/activityEstList(expIndex,2), activityEstList(expIndex,4)/activityEstList(expIndex,2), activityEstList(expIndex,6)/activityEstList(expIndex,2);];
    relFlightActDuration(expIndex,:)= [activityEstList(expIndex,2)/activityEstList(expIndex,2), activityEstList(expIndex,4)/activityEstList(expIndex,2), activityEstList(expIndex,6)/activityEstList(expIndex,4);];
 
    xNames{expIndex}= strcat(dataset(expIndex).fileName(1:8),' - (',dataset(expIndex).type, ' - ', dataset(expIndex).gender, ')');
end

%Plot the estimate flight actiovity
plotDataFlag=true;
if plotDataFlag
    figure()
    bp= bar(relFlightActDuration);
    %l= cell(1,3);
    %l{1}='Prior CO2'; l{2}='With CO2'; l{3}= 'Post CO2';
    %legend(l, 'Location', 'northeast');
    set(gca,'XtickLabel', xNames(2:2:end));
    if contains(version, 'R2019')
        xtickangle(90)    
    end
    %title('Estimated flight activity compared to prior-CO2 activity');
    title('Estimated flight activity compared to previous activity');
end

%% =============
% Delete from dataset the experiments where was a lower activity during the
% CO2 release than in their initial acclimatation hour.
% * Remember that you must move the h5 file deprecated to its corresponding
% deprecated folder (to don't load it again next time)
expDeprecated= find(relFlightActDuration(:,2) <= 1);
%Delete from the dataset the experiment not valid
if any(expDeprecated)
    %Build output file path
    fileName= strcat('Mosquito_Project_Report','.xlsx');
    sheetName= 'Experiments_Deprecated';
    outputFile= strcat(outputPath, fileName);
    %Check if the file exist
    if isfile(strcat(outputFile))
       %Load the name of the sheets existing on the file
       [t,sheetList]= xlsfinfo(outputFile);
       if any(strcmp(sheetList, sheetName))
           % if the sheet already exist in the file
           dataFromExcel = readtable(outputFile, 'sheet',sheetName);
           %Compare the exp. dates to see if the current exp. analysis are already on the Summary Report
           newEntries= setdiff(fileID, dataFromExcel.date);    
           if ~cellfun(@isempty,newEntries)
                % Find the indexes for newEntries experiments in fileNameList
                [newEntries,indexes1]= intersect(fileID, newEntries);               
                % Prepare the data to save in file: [expDate-expColorCueLabel-EstFlightAct_AIR-EstFlightAct_CO2-EstFlightAct_posCO2]
                expDeprecData=[fileID, outputFolder, num2cell(relFlightActDuration(indexes1,:))];
                writetable([dataFromExcel; expDeprecData], outputFile, 'sheet',sheetName);
            end
            clear dataFromExcel
       else
           % If there is no sheet called as sheetName
           sheetHeader= {'date', 'expColorCueLabel', 'EstFlightAct_AIR', 'EstFlightAct_CO2', 'EstFlightAct_posCO2'}
           expDeprecData=[fileNameList(expDeprecated,:), outputFolder, num2cell(relFlightActDuration(expDeprecated,:))];
           %add header to expDeprecData
           expDeprecData= [sheetHeader; expDeprecData];
           writecell(expDeprecData, outputFile, 'sheet',sheetName, 'Range', 'A1');
           % Clean temp variables from workspace
           clear sheetHeader expDeprecData
       end
    end
    % Move the file to its corresponding ExpDeprecated folder
    if ~isfolder(strcat(inputPath,'ExpDeprecated\'))
        mkdir(inputPath,'ExpDeprecated\');
    end
    movefile(filePath, strcat(inputPath,subFolder, 'ExpDeprecated\'));
    % Erase the experiment deprecated from the current dataset matrices
    clear all;
    %End code execution
    %return
    
    error('ERROR: The Experiment did not pass our control threshold (CO2 activity higher than initial AIR activity).');
    %dataset(expDeprecated)=[];
    %baseColorIndexList(expDeprecated)= [];
    %testedColorList(expDeprecated,:)= [];
    %filesList(expDeprecated)= []; 
    
end
% =============================================================

%% ===========================
% Plot the heatmaps for each of the experiments. There are 3 figures for
% each experiment heatmap prior/during/post CO2being released
% Plot the heatmap in for the trajectories inside the test section
if plotDataFlag
    %Find the indexes for AIR, CO2 and postCO2
    indexPrevCO2 = find(strcmp(dataset.stim(:), 'AIR'));
    indexWithCO2 = find(strcmp(dataset.stim(:),'CO2'));
    indexPostCO2 = find(strcmp(dataset.stim(:),'postCO2'));
    if dataset.type(1) =='m'
        insType= 'opsin_';
        if dataset.type(2) =='0'
            insType= strcat(insType,'1&2');    
        else
            insType= strcat(insType, dataset.type(2));
        end
    elseif dataset.type(1) == 'l'
        insType= strcat('line_',dataset.type(2));
    else
        insType= dataset.type;
    end   
    %Parameters related to the heatmaps
    topValueNormalized= 0.0001; %0.0003;
    imgTitle= strcat(dataset.fileName(1:15),'_',insType,'_',dataset.gender);
    plotCues= false;     %Control Flag to plot the odor/visual cues
    plot_XY_XZ_heatmaps_v5(plotCues,dataset.attr_x(indexPrevCO2), dataset.attr_y(indexPrevCO2), dataset.attr_z(indexPrevCO2), dataset.expCues, topValueNormalized, strcat('0_1_',imgTitle,'_heatmap_prevCO2'));    
    plot_XY_XZ_heatmaps_v5(plotCues,dataset.attr_x(indexWithCO2), dataset.attr_y(indexWithCO2), dataset.attr_z(indexWithCO2), dataset.expCues, topValueNormalized, strcat('0_2_',imgTitle,'_heatmap_withCO2'));
    plot_XY_XZ_heatmaps_v5(plotCues,dataset.attr_x(indexPostCO2), dataset.attr_y(indexPostCO2), dataset.attr_z(indexPostCO2), dataset.expCues, topValueNormalized, strcat('0_3_',imgTitle,'_heatmap_postCO2'));   
end
    % =============================================================


%% ===========================    
% Count the amount of trajectories detected while stim is prior, during and
% post CO2
trajSummary= zeros(length(dataset), 9);
%Find the indexes for AIR, CO2 and postCO2
indexPrevCO2 = find(strcmp(dataset.stim(:), 'AIR'));
indexWithCO2 = find(strcmp(dataset.stim(:),'CO2'));
indexPostCO2 = find(strcmp(dataset.stim(:),'postCO2'));
% Count the total number of trajectories longer than flightTimeLimit
% seconds. AS each row value in dataset.attr_id(:) is related to 1
% frame we only need to check in how many frames that ID has been
% detected
disp(strcat(' --- FILE: ',dataset.fileName, ' --- '));
disp(' --- Data prev CO2 --');  
[trajSummary(1,1), trajSummary(1,2), trajSummary(1,3)]= count_trajectories([dataset.attr_id(indexPrevCO2) dataset.attr_time(indexPrevCO2)], flightTimeLimit);
disp(' --- Data during CO2 --');
[trajSummary(1,4), trajSummary(1,5), trajSummary(1,6)]= count_trajectories([dataset.attr_id(indexWithCO2), dataset.attr_time(indexWithCO2)], flightTimeLimit);
disp(' --- Data post CO2 --');
[trajSummary(1,7), trajSummary(1,8), trajSummary(1,9)]= count_trajectories([dataset.attr_id(indexPostCO2), dataset.attr_time(indexPostCO2)], flightTimeLimit);
disp(' =====================================');
% =======================
% =======================

%% Generate your own speed calculator

for odorChecked=[{'AIR'},{'CO2'},{'postCO2'}]
    for expIndex=1:length(dataset)
        odorIndex= find(strcmp(dataset(expIndex).stim, odorChecked));
        data=[dataset(expIndex).attr_id(odorIndex), dataset(expIndex).attr_time(odorIndex), dataset(expIndex).attr_x(odorIndex), dataset(expIndex).attr_y(odorIndex), dataset(expIndex).attr_z(odorIndex)];
        %select all the IDs from the vector
        uniqueID=unique(data(:,1));
        %Transpose from a column matrix to a row matrix
        uniqueID=uniqueID';
        % Initialice local variables
        totalTraj=0; 
        flightTimeList= nan(length(uniqueID),1);
        totalFlightTime=0;
        speedList=[];
        objIDsList=[];
        % For each insect ID, check if 
        for index= 1:length(uniqueID)
            objID= uniqueID(index);
            %Load the frames where appears the current objID
            objTime= find(data(:,1) == objID);
            % Estimate the duration of each insect ID trajectory
            duration= get_trajectory_duration(data(objTime(:),2));

            %Consider only if flight duration is bigger than flightTimeLimit seconds
            if duration >= flightTimeLimit
                % Calculate the lenght traveled by trajectory ID in each axis
                distX= sum(abs(diff(data(objTime(:),3))));
                distY= sum(abs(diff(data(objTime(:),4))));
                distZ= sum(abs(diff(data(objTime(:),5))));

                % Estimate the mean velocity
                distT= sqrt(distX^2 + distY^2 + distZ^2);
                speed= distT/duration;

                speedList= [speedList, speed]; 
                objIDsList= [objIDsList, objID];
            end
        end
        
        disp(strcat('-- Iteration: ', int2str(expIndex), ' ---- ', odorChecked, ' -- max speed:', sprintf('%.6f',max(speedList))))
        
        % Store the mean speed of each trajectory ID
        if strcmp(odorChecked, 'AIR')
            meanSpeedList(expIndex, 1)= mean(speedList);
            meanSpeedExp(expIndex).IDsAIR= objIDsList;
            meanSpeedExp(expIndex).speedsAIR= speedList;
        elseif strcmp(odorChecked, 'CO2')
            meanSpeedList(expIndex, 2)= mean(speedList);
            meanSpeedExp(expIndex).IDsCO2= objIDsList;
            meanSpeedExp(expIndex).speedsCO2= speedList;
        elseif strcmp(odorChecked, 'postCO2')
            meanSpeedList(expIndex, 3)= mean(speedList);
            meanSpeedExp(expIndex).IDsPostCO2= objIDsList;
            meanSpeedExp(expIndex).speedsPostCO2= speedList;
        end       

    end
end 

plot(meanSpeedList, '*')
lg=[{'AIR - CO2 - postCO2'}]
legend(lg, 'Location', 'northeast');

%% ==========================
% Choose a random spot in the wind tunnel to compare its number of visits
% with the base and test cues volumes in the next section below this one.
% The same random volume will be used for AIR and CO2 analysis for the same
% experiment

% Create a random position for a 3rd volume inside the wind tunnel
radius= 0.07; % Volume radius
%Pick a random spot inside the XY plane of the test section (the volume
%must be fully inside the test section
cuesX= -lim_x + cell2mat(dataset.expCues(2,2));
tmpLimX= lim_x - radius;
tmpLimY= lim_y - radius;
tmpX = (cuesX+radius) + rand*(tmpLimX-(cuesX+radius));
tmpY = -tmpLimY + rand*(tmpLimY-(-tmpLimY));

centerList(1,:)= [tmpX tmpY];   


%% ===========================
% % => (A.) This loop runs for each type of mosquito (wt,mX, lX), gender (m, f) or odor used (AIR, CO2, postCO2) <= 
% % => To generate several datafiles, run it several times changing the values checkedType, checkedGender and checkedOdor <=  
checkedOdor= 'CO2';
%dataFolder= 'analysisData\';
dataFolder= 'analysisData_newVersion\';
xNames= cell(1,length(dataset));
countListWithOdor= zeros(1,length(dataset.expCues(2:end,1))+1);


% disp(strcat(' - Working with expIndex:',num2str(expIndex)));
% Create the new file Excel file name
disp(strcat(' * Working with group: ',{' '}, dataset.fileName(1:15)));
fileName= strcat(dataset.fileName(1:15), '_countsInsideCueVol_',checkedOdor,'.xlsx');
% Returns vector with the blobal counts per each cue position [counts in pos1, counts in pos2, ..., counts in posN]
if strcmp(dataFolder,'analysisData\') 
    countListWithOdor(1,:)= count_insect_in_volume_v6(strcat(outputPath, outputFolder, dataFolder), fileName, dataset, flightTimeLimit, checkedOdor);
else
    disp('Using new version')
   % New version. It load the subset related to the odorChecked and then do the counts
   % Load Indexes
   indexWithOdor=[];
   indexWithOdor = find(strcmp(dataset.stim(:),checkedOdor));
   dataEntry.expCues= dataset.expCues;
   dataEntry.attr_id= dataset.attr_id(indexWithOdor);
   dataEntry.attr_time= dataset.attr_time(indexWithOdor);
   dataEntry.attr_x= dataset.attr_x(indexWithOdor);
   dataEntry.attr_y= dataset.attr_y(indexWithOdor);
   dataEntry.attr_z= dataset.attr_z(indexWithOdor);
   % Counts elements near the cues for the subset related to the odorChecked 
   countListWithOdor(1,:)= count_insect_in_volume_v7(strcat(outputPath, outputFolder, dataFolder), fileName, dataEntry, flightTimeLimit, checkedOdor, centerList(1,:));
end
% Pick the experiment dates to plot later as the X axis legend
%xNames{expIndex}= strcat(dataset(expIndex).fileName(1:8),' - (',dataset(expIndex).type, ' - ', dataset(expIndex).gender, ')');
clear checkedType checkedGender checkedOdor dataFolder
% => END (A.) LOOP <=




%% ===========================
% PI FOR 2 COLORS ONLY!
% Generates the PreferenceIndex value per each experiment date

piList= generate_PI(countListWithOdor(1:2), baseColorIndex);
%Change the NaN by 0
piList(isnan(piList))=0;

if plotDataFlag
    figure()
    bar(piList)
    ylim([-1 1]);
    set(gca,'XtickLabel', xNames);
    
    if contains(version, 'R2019')
        xtickangle(45)    
    end
    title('Preference Index towards the tested color');
end


%% =============== ADD ROW TO REPORT ===============
% Add a row with the iformation from the analysis into the report
baseColor= 'white';  %use for exp red_vs_white
fileName= strcat('Mosquito_Project_Report','.xlsx');
%sheetName= strcat('exp_with_base_color_',baseColor);
sheetName= outputFolder(1:end-1);

sheetName= strcat('MUTANTS_with_base_color_',baseColor);
fileName= strcat('JeffTest_Data_Mosquito_Mutants_Project_Report','.xls');

% Add the total duration of all trajectories detected with an odor
% (activityEstList) in the trajSummary variable.
loadTotalFlighDurationPerODor= false;
if loadTotalFlighDurationPerODor
    % If the report contain a column with the sum of all trajectories ID
    % duration per odor
    fullTrajSummary= [activityEstList(:,1:2),trajSummary(:,2:3),activityEstList(:,3:4),trajSummary(:,5:6),activityEstList(:,5),trajSummary(:,8:9)];
else
    fullTrajSummary= trajSummary;
end    
add_entry_in_summary_report(outputPath, fileName, sheetName, fileID, fullTrajSummary, countListWithOdor(:,1:2), baseColorIndex, testedColor, piList)





%% ===========================
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

