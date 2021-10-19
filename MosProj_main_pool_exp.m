
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

%% (0.) ========== USER PARAMETERS ==========
close all;
clear all;
format long;
%opengl('save', 'software');
% Information related to the data to analyze
%workspace = 'C:\Users\Riffelllab\Documents\Matlab\AnalysisFLYDRA\';
%workspace= 'C:\Users\dalon\Documents\MATLAB\AnalysisFLYDRA_v011920\';
workspace= 'C:\Users\dalon\Documents\GitHub\AnalysisFLYDRA_2020';
inputPath= 'D:\DataFlydra\MosquitoProject\INPUT_Data\';
outputPath= 'D:\DataFlydra\MosquitoProject\OUTPUT_Graphs\';
% Variables to use in this script
fps=60.0; %90.0;
% Lower time threshold (in seconds) to plot a trajectory or not
flightTimeLimit= 1.5;
%Select the Control Color value
baseColor= 'white'; %select between: white, R-HUE, GcT1
% Test section X,Y,Z limits (start_x= -lim_x and end_x= +lim_x --> Total x distance is lim_x*2)
% These lim values are the current size of the Test Section (DO NOT CHANGE IT)
lim_x= 0.9144;
lim_y= 0.3048;
lim_z= 0.6096;

%subFolder= 'FLYDRA_Trajectories_4_blacks_squared\';
%subFolder= 'FLYDRA_Trajectories_4_colors_squared\';
%subFolder= 'FLYDRA_Trajectories_4_reds_squared\';
%subFolder= 'FLYDRA_Trajectories_4_grays_squared\';

%subFolder= 'FLYDRA_Trajectories_gray_red\';
%subFolder= 'FLYDRA_Trajectories_black_red\';


subFolder= 'FLYDRA_Trajectories_white_black\';
%subFolder= 'FLYDRA_Trajectories_white_blue\';
%subFolder= 'FLYDRA_Trajectories_white_gray\';
%subFolder= 'FLYDRA_Trajectories_white_green\';
%subFolder= 'FLYDRA_Trajectories_white_purple\';
%subFolder= 'FLYDRA_Trajectories_white_red\';
%subFolder= 'FLYDRA_Trajectories_white_yellow\';
%subFolder= 'FLYDRA_Trajectories_white_pR01\';

%subFolder= 'FLYDRA_Trajectories_2_blacks\';
%subFolder= 'FLYDRA_Trajectories_gct1_gwt1\';
%subFolder= 'FLYDRA_Trajectories_mutants_lineX\';

%subFolder= 'FLYDRA_Trajectories_AnStephensi\';
%subFolder= 'FLYDRA_Trajectories_CxQuinquefasciatus\';
%subFolder= 'FLYDRA_Trajectories_HistamineFed\';

%Create the outputFolder folder
outputFolder = load_output_folder(outputPath, subFolder);
% ===========================
% ===========================



%% ------------------------------------------------------------------------
% METHODS TO GENERATE SUMMARIES/FILES FROM RAW DATA
% -------------------------------------------------------------------------

%% (1.) Load the name of the files containing the experiments datasets stored in subFolder
inputPath= strcat(inputPath,subFolder);
cd(inputPath)
filesList=dir('*mainbrain.h5');

% Move to the MATLAB workspace
cd(workspace);

%initialize list for tested color used and position ob tyhe base color cue in WT
baseColorIndexList=zeros(length(filesList),1);
testedColorList= cell(length(filesList),1);

loadFullDataset=true; % Variable used to test with partial data
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
        indexWithCO2= find(dataset(expIndex).attr_time(:) >= ts_startCO2 & dataset(expIndex).attr_time(:) < ts_endCO2);
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
clear fileName loadFullDataset indexPrevCO2 indexWithCO2 indexPostCO2
clear attr_id attr_time attr_frame attr_x attr_y attr_z
clear cuesSetup ts_startAIR ts_startCO2 ts_endCO2 ts_endAIR mType mGender
clear filesList expIndexes 
% ===========================
% ===========================


%% (2) ============ ESTIMATE FLIGHT ACTIVITY ===============
% Estimate relative activity by adding all trajectories duration and then
% dividing it by the total number of trajectories. (different estimations 
% for prior/during/post CO2) 
% WARNING: WE are dividing each of the total number of trajectories prior,
% during and post by the total number of trajectories prior CO2
activityEstList= zeros(length(dataset), 6);
relFlightActDuration= zeros(length(dataset), 3);
xNames= cell(1,length(dataset));
for expIndex=1:length(dataset)
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

figure()
bar(relFlightActDuration(1:end,:));
l= cell(1,3);
l{1}='Prior CO2'; l{2}='With CO2'; l{3}= 'Post CO2';
legend(l, 'Location', 'northeast');
if length(dataset) > 10
    set(gca,'XtickLabel', xNames(2:2:end));
else
    set(gca,'XtickLabel', xNames);
end
if contains(version, 'R2019')
    xtickangle(90)    
end
%title('Estimated flight activity compared to prior-CO2 activity');
title('Estimated flight activity compared to previous activity');
% ===========================
% ===========================


%% (3.)========== DELETE EXP. THAT FAIL CONTROL TRESHOLD =============
% Delete from dataset the experiments where was a lower activity during the
% CO2 release than in their initial acclimatation hour.
% * Remember that you must move the h5 file deprecated to its corresponding
% deprecated folder (to don't load it again next time)
indexExpDeprecated= find(relFlightActDuration(:,2) <= 1);
%Delete from the dataset the experiment not valid
if any(indexExpDeprecated)
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
           newEntries= setdiff(fileNameList(indexExpDeprecated), dataFromExcel.date);    
           if cellfun(@any,newEntries)
                % Find the indexes for newEntries experiments in fileNameList
                %[newEntries,indexes1]= intersect(fileNameList(indexExpDeprecated), newEntries);               
                % Prepare the data to save in file: [expDate-expColorCueLabel-EstFlightAct_AIR-EstFlightAct_CO2-EstFlightAct_posCO2
                [tempOutFolderList{1:length(newEntries),1}]= deal(outputFolder);
                expDeprecData=[fileNameList(indexExpDeprecated,:), tempOutFolderList, num2cell(relFlightActDuration(indexExpDeprecated,:))];
                writetable([dataFromExcel; expDeprecData], outputFile, 'sheet',sheetName);
            end
            clear dataFromExcel newEntries
       else
           % If there is no sheet called as sheetName
           sheetHeader= {'date', 'expColorCueLabel', 'EstFlightAct_AIR', 'EstFlightAct_CO2', 'EstFlightAct_posCO2'};
           expDeprecData=[fileNameList(indexExpDeprecated,:), outputFolder, num2cell(relFlightActDuration(indexExpDeprecated,:))];
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
    for i= 1:length(indexExpDeprecated)
        movefile(strcat(inputPath,dataset(indexExpDeprecated(i)).fileName), strcat(inputPath,'ExpDeprecated\'));
    end
    % Erase the experiment deprecated from the current dataset matrices
    dataset(indexExpDeprecated)=[];
    baseColorIndexList(indexExpDeprecated)= [];
    testedColorList(indexExpDeprecated,:)= [];
    fileNameList(indexExpDeprecated)= []; 
    relFlightActDuration(indexExpDeprecated,:)= [];
    activityEstList(indexExpDeprecated, :)=[];
    xNames(indexExpDeprecated)=[];
end
% ===========================
% ===========================


%% (4.) ========== COUNT TRAJECTORIES DETECTED  =================    
% Count the amount of trajectories detected while stim is prior, during and
% post CO2
trajSummary= zeros(length(dataset), 9);
for expIndex= 1:length(dataset)
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
clear indexPrevCO2 indexWithCO2 indexPostCO2
% ===========================
% ===========================

%% (4.1) ========== CALCULATE MEAN TRAJECTORIES SPEED  =================    

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

            end
        end
        
        disp(strcat('-- ', int2str(expIndex), ' ---- ', odorChecked))
        disp(strcat('max speed:', sprintf('%.6f',max(speedList))))
        if strcmp(odorChecked, 'AIR')
            meanSpeedList(expIndex, 1)= mean(speedList);
        elseif strcmp(odorChecked, 'CO2')
            meanSpeedList(expIndex, 2)= mean(speedList);
        elseif strcmp(odorChecked, 'postCO2')
            meanSpeedList(expIndex, 3)= mean(speedList);
        end
    end
end 

plot(meanSpeedList)
lg=[{'AIR'}, {'CO2'}, {'postCO2'}]
legend(lg, 'Location', 'northeast');
% ===========================
% ===========================



%% (5.) ======== PLOT/SAVE HEATMAPS ===================
% Plot the heatmaps for each of the experiments. There are 3 figures for
% each experiment heatmap prior/during/post CO2being released
% Heatmpas control parameters
topValueNormalized= 0.0001; %0.0002;
plotCues= false;     %Control Flag to plot the odor/visual cues
    
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
        insType= dataset(expIndex).type;
    end
    imgTitle= strcat(dataset(expIndex).fileName(1:15),'_',insType,'_',dataset(expIndex).gender);
    plot_XY_XZ_heatmaps_v5(plotCues,dataset(expIndex).attr_x(indexPrevCO2), dataset(expIndex).attr_y(indexPrevCO2), dataset(expIndex).attr_z(indexPrevCO2), dataset(expIndex).expCues, topValueNormalized, strcat(num2str(expIndex),'_1_',imgTitle,'_heatmap_prevCO2'));
    plot_XY_XZ_heatmaps_v5(plotCues,dataset(expIndex).attr_x(indexWithCO2), dataset(expIndex).attr_y(indexWithCO2), dataset(expIndex).attr_z(indexWithCO2), dataset(expIndex).expCues, topValueNormalized, strcat(num2str(expIndex),'_2_',imgTitle,'_heatmap_withCO2'));
    plot_XY_XZ_heatmaps_v5(plotCues,dataset(expIndex).attr_x(indexPostCO2), dataset(expIndex).attr_y(indexPostCO2), dataset(expIndex).attr_z(indexPostCO2), dataset(expIndex).expCues, topValueNormalized, strcat(num2str(expIndex),'_3_',imgTitle,'_heatmap_postCO2'));   
end    
clear indexPrevCO2 indexWithCO2 indexPostCO2 imgTitle insType plotCues
% ===========================
% ===========================

%% (6.) ======== TRAJECTORIES SPEED CALCULATOR ===================
%  Estimate the mean trajectories speed by dividing the distance traveled
%  by the trajectory duration

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

            end
        end
        
        disp(strcat('-- ', int2str(expIndex), ' ---- ', odorChecked))
        disp(strcat('max speed:', sprintf('%.6f',max(speedList))))
        if strcmp(odorChecked, 'AIR')
            meanSpeedList(expIndex, 1)= mean(speedList);
        elseif strcmp(odorChecked, 'CO2')
            meanSpeedList(expIndex, 2)= mean(speedList);
        elseif strcmp(odorChecked, 'postCO2')
            meanSpeedList(expIndex, 3)= mean(speedList);
        end
    end
end 

plot(meanSpeedList)
lg=[{'AIR'}, {'CO2'}, {'postCO2'}]
legend(lg, 'Location', 'northeast');




%% (7.) ========= CREATE RANDOM VOLUMES POSITIONS ==================
% Choose a random spot in the wind tunnel to compare its number of visits
% with the base and test cues volumes in the next section below this one.
% The same random volume will be used for AIR and CO2 analysis for the same
% experiment

% Create a random position for a 3rd volume inside the wind tunnel
radius= 0.07; % Volume radius

for expIndex=1: length(dataset)
    %Pick a random spot inside the XY plane of the test section (the volume
    %must be fully inside the test section
    cuesX= -lim_x + cell2mat(dataset(expIndex).expCues(2,2));
    tmpLimX= lim_x - radius;
    tmpLimY= lim_y - radius;
    tmpX = rand*(tmpLimX);
    tmpY = -tmpLimY + rand*(tmpLimY-(-tmpLimY));

    centerList(expIndex,:)= [tmpX tmpY];   
end
% ===========================
% ===========================



%% (8.) ========== COUNT VISIT NEAR VOLUMES (TEST/BASE/RANDOM) =================
% Count insects (of particular type and gender) near visualClue when the CO2 is being released
% An insect ID will be counted as many times as it passes over the volme
% 
% % => (7.A.) This loop runs for each type of mosquito (wt,mX, lX), gender (m, f) or odor used (AIR, CO2, postCO2) <= 
% % => To generate several datafiles, run it several times changing the values checkedType, checkedGender and checkedOdor <=  
checkedOdor= 'CO2';
onlyOneHourCO2= false;   %To use in experiments with 1 and 2 hours of CO2
%dataFolder= 'analysisData_newVersion_1stHourOnly\';

%dataFolder= 'analysisData\';
dataFolder= 'analysisData_newVersion\';
xNames= cell(1,length(dataset));
% countListWithOdor has info from cues + 1 random volume (centerList(:)
countListWithOdor= zeros(length(dataset),length(dataset(1).expCues(2:end,1))+1);


for expIndex= 1:length(dataset)
    % disp(strcat(' - Working with expIndex:',num2str(expIndex)));
    % Create the new file Excel file name
    disp(strcat(' * Working with group: ',{' '}, dataset(expIndex).fileName(1:15)));
    fileName= strcat(dataset(expIndex).fileName(1:15), '_countsInsideCueVol_',checkedOdor,'.xlsx');
    % Returns vector with the blobal counts per each cue position [counts in pos1, counts in pos2, ..., counts in posN]
    if strcmp(dataFolder,'analysisData\') 
        countListWithOdor(expIndex,:)= count_insect_in_volume_v6(strcat(outputPath, outputFolder, dataFolder), fileName, dataset(expIndex), flightTimeLimit, checkedOdor);
    else
        disp('Using new version')
       % New version. It load the subset related to the odorChecked and then do the counts
       % Load Indexes
       indexWithOdor=[];
       indexWithOdor = find(strcmp(dataset(expIndex).stim(:),checkedOdor));
       if onlyOneHourCO2 && (dataset(expIndex).attr_time(indexWithOdor(end)) - dataset(expIndex).attr_time(indexWithOdor(1)) > 3600)
           disp(' - Analyzing only for the first hour of CO2');
           endFstHourOdor= find(dataset(expIndex).attr_time(indexWithOdor) > dataset(expIndex).attr_time(indexWithOdor(1))+3600,1);
           indexWithOdor= indexWithOdor(1:endFstHourOdor); 
       end
           
       dataEntry.expCues= dataset(expIndex).expCues;
       dataEntry.attr_id= dataset(expIndex).attr_id(indexWithOdor);
       dataEntry.attr_time= dataset(expIndex).attr_time(indexWithOdor);
       dataEntry.attr_x= dataset(expIndex).attr_x(indexWithOdor);
       dataEntry.attr_y= dataset(expIndex).attr_y(indexWithOdor);
       dataEntry.attr_z= dataset(expIndex).attr_z(indexWithOdor);
       % Counts elements near the cues for the subset related to the odorChecked 
       countListWithOdor(expIndex,:)= count_insect_in_volume_v7(strcat(outputPath, outputFolder, dataFolder), fileName, dataEntry, flightTimeLimit, checkedOdor, centerList(expIndex,:));
    end
end
clear checkedType checkedGender checkedOdor dataFolder dataEntry
% => END (7.A.) LOOP <=
% ===========================
% ===========================


%% (9.) ========GENERATE PI FOR 2 COLORS ONLY (TEST vs BASE) ==================
% From the # of visits generated in 7.A Section, Generate the PreferenceIndex value per each experiment date
piList= generate_PI(countListWithOdor(:,1:2), baseColorIndexList);
%Change the NaN by 0
piList(isnan(piList))=0;
figure()
bar(piList)
ylim([-1 1]);
set(gca,'XtickLabel', xNames);

if contains(version, 'R2019')
    xtickangle(45)    
end
title('Preference Index towards the tested color');
% ===========================
% ===========================



%% (10.) =============== ADD DATA IN SUMMARY REPORT ===============
% generate report with all the information for this type of experiments
%countListWithOdor= zeros(length(dataset),length(dataset(1).expCues(2:end,1)));
%baseCueIndexList=zeros(length(dataset(:)),1);
%baseColor='None';
%piList=0;
%baseColor= 'R-HUE';  %use for exp red_vs_white
fileName= 'Mosquito_Project_Report.xlsx';
%sheetName= 'AnStephensi_exp';
%sheetName= 'CxQuinquefasciatus_exp';
% for color1 vs color2 experiments (no common base color between experiment groups)
%sheetName= outputFolder(1:end-1);
sheetName= strcat('exp_with_base_color_',baseColor);

% Add the total duration of all trajectories detected with an odor
% (activityEstList) in the trajSummary variable.
fullTrajSummary= [activityEstList(:,1:2),trajSummary(:,2:3),activityEstList(:,3:4),trajSummary(:,5:6),activityEstList(:,5:6),trajSummary(:,8:9)];

if length(dataset) == 1 
    %From the group of experiments there is only one that passed our initial activity control
    add_entry_in_summary_report(outputPath, fileName, sheetName, fileNameList, fullTrajSummary, countListWithOdor(:,1:2), baseColorIndexList, testedColorList, piList)
else
    % From the group of experiments, there are several that passed our initial activity control
    generate_report(outputPath, fileName, sheetName, fileNameList, fullTrajSummary, countListWithOdor(:,1:2), baseColorIndexList, testedColorList, piList)
end
% ===========================
% ===========================




%% ------------------------------------------------------------------------
% METHODS TO WORK WITH ALL THE DATA PREVIOUSLY GENERATED IN SUMMARIES/FILES
% -------------------------------------------------------------------------

%% (0.)========== USER PARAMETERS ==========
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
% Variables to use in this script
fps=60.0; %90.0;
% Lower time threshold (in seconds) to plot a trajectory or not
flightTimeLimit= 1.5;
%Select the Control Color value
baseColor= 'white'; %select between: white, R-HUE, GcT1
% Test section X,Y,Z limits (start_x= -lim_x and end_x= +lim_x --> Total x distance is lim_x*2)
% These lim values are the current size of the Test Section (DO NOT CHANGE IT)
lim_x= 0.9144;
lim_y= 0.3048;
lim_z= 0.6096;

%subFolder= 'FLYDRA_Trajectories_4_blacks_squared\';
%subFolder= 'FLYDRA_Trajectories_4_colors_squared\';
%subFolder= 'FLYDRA_Trajectories_4_reds_squared\';
%subFolder= 'FLYDRA_Trajectories_4_grays_squared\';

%subFolder= 'FLYDRA_Trajectories_gray_red\';
%subFolder= 'FLYDRA_Trajectories_black_red\';


%subFolder= 'FLYDRA_Trajectories_white_black\';
%subFolder= 'FLYDRA_Trajectories_white_blue\';
%subFolder= 'FLYDRA_Trajectories_white_gray\';
subFolder= 'FLYDRA_Trajectories_white_green\';
%subFolder= 'FLYDRA_Trajectories_white_purple\';
%subFolder= 'FLYDRA_Trajectories_white_red\';
%subFolder= 'FLYDRA_Trajectories_white_yellow\';
%subFolder= 'FLYDRA_Trajectories_white_pR01\';



%subFolder= 'FLYDRA_Trajectories_2_blacks\';
%subFolder= 'FLYDRA_Trajectories_gct1_gwt1\';
%subFolder= 'FLYDRA_Trajectories_mutants_lineX\';

%subFolder= 'FLYDRA_Trajectories_AnStephensi\';
%subFolder= 'FLYDRA_Trajectories_CxQuinquefasciatus\';
%subFolder= 'FLYDRA_Trajectories_HistamineFed\';

%Create the outputFolder folder
outputFolder = load_output_folder(outputPath, subFolder);
% ===========================
% ===========================


%% (11.) =============== LOAD SUMMARY REPORT ===============
% Load the report with all the information for all experiments with a sheet
% in the file. The function returns the following tables:
% - smryTrajectories= [date numOfTrajAIR avgTimeFltAIR maxTimeFltAIR ..
%                       numOfTrajCO2 avgTimeFltCO2 maxTimeFltCO2 numOfTrajPostCO2 avgTimeFlt_PostCO2 maxTimeFltPostCO2]
% - smryCounts=       [date testedColor Pos-1-farthest Pos-2-closest baseColor-Pos]
% - smryPIs=          [date PI towards color]

% Choose the name of the xlsx sheet to load the data into the workspace
sheetName= strcat('exp_with_base_color_', baseColor);
%sheetName= outputFolder(1:end-1);
[smryTrajectories, smryCounts, smryPIs]= load_summary_report(outputPath, 'Mosquito_Project_Report.xlsx', sheetName);
clear sheet_name
% ===========================
% ===========================



%% (12.) =================== WORKING WITH INDIVIDUAL FILES PER EXPERIMENT ===================
%  - CALCULATE PI VALUES FOR EACH TRAJ ID INSIDE VOLUME (CONTROL/TEST/RANDOM)
% From the summary report file, find the indexes for the experiments we
% want to analyze ('gray2.5', 'gray4.0', 'gray4.5', 'gray6.5', 'gray9.5')
%(1.) Estimate time spent near visual cues for each one of the trajectories

%Select the the odor and if working with all experiments or a particular
%color subset
odorChecked='CO2';
workWithAllColors= true;

% Select with which data do you want to work (1 color/group or all colors) 
if workWithAllColors
    baseColorIndexList= smryCounts.baseColor_Pos;
    testColorList= unique(smryCounts.testedColor);
    % Load all files from group folder
    %filesPath= strcat(outputPath, 'AnalysisData_All_Exp_BaseC_Red\');
    %filesPath= strcat(outputPath, 'AnalysisData_WithRndVol_All_Exp_baseC_rhue\');
    filesPath= strcat(outputPath, 'AnalysisData_ForPaper_All_Exp_BaseC_White\');
    % Load files names
    cd(filesPath);
    tmp=dir(strcat('*countsInsideCueVol_',odorChecked,'.xlsx'));
    filesList= vertcat(tmp.name);
    clear tmp
else
    if contains(outputFolder, 'AnStephensi') || contains(outputFolder, 'CxQuinquefasciatus')
        baseColorIndexList= smryCounts.baseColor_Pos;
        testColorList= unique(smryCounts.testedColor);
        % Load the files in the individual folder
        filesPath= strcat(outputPath, outputFolder,'analysisData_newVersion\');
        % Load files names
        cd(filesPath);
        tmp=dir(strcat('*countsInsideCueVol_',odorChecked,'.xlsx'));
        filesList= vertcat(tmp.name);
        clear tmp
        
    else
        colorGroup= 'black';   %Select the color to work with
        cIndexes= find(strncmp(smryCounts.testedColor, colorGroup, 4));
        baseColorIndexList= smryCounts.baseColor_Pos(cIndexes);
        testColorList= unique(smryCounts.testedColor(cIndexes));
        % Load the files in the individual folder
        filesPath= strcat(outputPath, outputFolder,'analysisData_newVersion\');
        %filesPath= strcat(outputPath, outputFolder,'analysisData_newVersion_h016\');
        % Load files
        cd(filesPath);
        filesList= strcat(char(smryCounts.date(cIndexes)),'_countsInsideCueVol_',odorChecked,'.xlsx');
       
    end
end
   

% Load the data from the set of files we are working
cd(workspace);
for fileIndex= 1:size(filesList,1)
    fileName= filesList(fileIndex,:);
    disp(strcat(' - Working with file: ', {' '}, fileName));
    filesName(fileIndex)= {fileName(1:8)};
    % data loaded structure: [vCuePosition	objID	timeStamp NaN initialTimeStamp finalTimeStamp]
    dataFromExcel = table2array(readtable(strcat(filesPath, fileName)));
    expData(fileIndex).name= fileName(1:15);
    % Initialize the correct fields in function of the odor used
    if strcmp(odorChecked, 'AIR')
        expData(fileIndex).listTotalIDsNearCueAIR=[];
        expData(fileIndex).PrefIndexTotalIDsAIR=[]; 
        expData(fileIndex).listTotalIDsRndVolAIR=[];
        expData(fileIndex).PrefIndexTotalIDsRndVolAIR=[];       
    elseif strcmp(odorChecked, 'CO2')
        expData(fileIndex).listTotalIDsNearCueCO2=[];
        expData(fileIndex).PrefIndexTotalIDsCO2=[];
        expData(fileIndex).listTotalIDsRndVolCO2=[];
        expData(fileIndex).PrefIndexTotalIDsRndVolCO2=[];   
    elseif strcmp(odorChecked, 'postCO2')
        expData(fileIndex).listTotalIDsNearCuePostCO2=[];
        expData(fileIndex).PrefIndexTotalIDsPostCO2=[];
        expData(fileIndex).listTotalIDsRndVolPostCO2=[];
        expData(fileIndex).PrefIndexTotalIDsRndVolPostCO2=[];   
    
    end
    if isempty(dataFromExcel)
        disp(strcat('   -> No counts detected in file: ',{' '}, fileName));
    else
        %Initialize the object that will contain the information for the current experiment
        tmpPIs=[];
        tmpIDs=[];
        tmpIDs1=[];
        tmpIDs2=[];
        tmpIDs3=[];
        tmpIDsInBoth=[];
        tmpPIsInRnd=[];
        tmpIDsInBaseAndRnd=[]; 
        % To load timestamp data (time spent by insect near the cue
        % Sorted the the data by cue visited ID timestamp
        sortedData= sortrows(dataFromExcel(:,1:3));
        if length(unique(sortedData(:,1)))== 1
            % Find all the Traj IDs in the file
            tmpIDs= unique(sortedData(:,2));
            if baseColorIndexList(fileIndex) == sortedData(1,1)
                % All the counts where detected near the Base color cue (PI == -1 for all traj IDs)
                tmpPIs= ones(length(tmpIDs),1)*(-1);
            else
                % All the counts where detected near the Tested color cue (PI == 1 for all traj IDs)
                tmpPIs= ones(length(tmpIDs),1);
            end
            %Add the IDs and the their time spent visiting the cue to the structure
            if strcmp(odorChecked, 'AIR')
                expData(fileIndex).listTotalIDsNearCueAIR=tmpIDs';
                expData(fileIndex).PrefIndexTotalIDsAIR=tmpPIs;
            elseif strcmp(odorChecked, 'CO2')           
                expData(fileIndex).listTotalIDsNearCueCO2=tmpIDs';
                expData(fileIndex).PrefIndexTotalIDsCO2=tmpPIs';
            elseif strcmp(odorChecked, 'postCO2')           
                expData(fileIndex).listTotalIDsNearCuePostCO2=tmpIDs';
                expData(fileIndex).PrefIndexTotalIDsPostCO2=tmpPIs';
            
            end
        else  
            %Group the IDs visitin each volume and count the time spent
            indexesPos= find(sortedData(:,1)==1);
            % Count # of frames where each Traj ID visit cue 1
            [tmpCts1, tmpIDs1]= groupcounts(sortedData(indexesPos,2));
            indexesPos= find(sortedData(:,1)==2);
            % Count # of frames where each Traj ID visit cue 2
            [tmpCts2, tmpIDs2]= groupcounts(sortedData(indexesPos,2));        
            indexesPos= find(sortedData(:,1)==3);
            % Count # of frames where each Traj ID visit cue 2
            [tmpCts3, tmpIDs3]= groupcounts(sortedData(indexesPos,2));  

            % check if any ID visited Base & Test or Base & random cues/vol
            if baseColorIndexList(fileIndex)== 1         
                [tmpIDsInBoth, inTmpIDs1, inTmpIDs2]= intersect(tmpIDs1, tmpIDs2);
            else
                [tmpIDsInBoth, inTmpIDs1, inTmpIDs2]= intersect(tmpIDs1, tmpIDs2);
            end
            
            % Generate the PI values for Tested Cue vs Base Cue
            if any(tmpIDsInBoth)
                % Assign PI values in function of where the Base cue was placed                
                if baseColorIndexList(fileIndex)== 1
                    %Base color was placed in position 1 (+Y axis) in current experiment
                    for i= 1:length(tmpIDsInBoth)
                        tmpPIs(i,1)= (tmpCts2(inTmpIDs2(i))-tmpCts1(inTmpIDs1(i)))/(tmpCts2(inTmpIDs2(i))+tmpCts1(inTmpIDs1(i)));
                    end
                    % Remove ID visiting both cues as they are grouped in tmpIDsInBoth
                    % (Erase duplicate ID when it has visited both cues)
                    tmpIDs1= setdiff(tmpIDs1, tmpIDsInBoth);
                    tmpIDs2= setdiff(tmpIDs2, tmpIDsInBoth);
                    % Assign their Pref Index value
                    tmpPIsOnly1= ones(length(tmpIDs1),1)*(-1);  %PI= -1
                    tmpPIsOnly2= ones(length(tmpIDs2),1);       %PI= 1
                else
                    %Base color was placed in position 2 (-Y axis) in current experiment
                    for i= 1:length(tmpIDsInBoth)
                        tmpPIs(i,1)= (tmpCts1(inTmpIDs1(i))-tmpCts2(inTmpIDs2(i)))/(tmpCts1(inTmpIDs1(i))+tmpCts2(inTmpIDs2(i)));
                    end
                    % Remove ID visiting both cues as they are grouped in tmpIDsInBoth
                    % (Erase duplicate ID when it has visited both cues)
                    tmpIDs2= setdiff(tmpIDs2, tmpIDsInBoth);
                    tmpIDs1= setdiff(tmpIDs1, tmpIDsInBoth);
                    tmpPIsOnly1= ones(length(tmpIDs1),1);       %PI= 1
                    tmpPIsOnly2= ones(length(tmpIDs2),1)*(-1);  %PI= -1                                       
                end
            else
                % ANY trajectory have visited BOTH cues
                if baseColorIndexList(fileIndex)== 1
                    %Base color was placed in position 1 (+Y axis) in current experiment
                    tmpPIsOnly1= ones(length(tmpIDs1),1)*(-1);  %PI= -1
                    tmpPIsOnly2= ones(length(tmpIDs2),1);       %PI= 1
                else
                    %Base color was placed in position 2 (-Y axis) in current experiment
                    tmpPIsOnly1= ones(length(tmpIDs1),1);       %PI= 1
                    tmpPIsOnly2= ones(length(tmpIDs2),1)*(-1);  %PI= -1
                end
            end
            
            %Add the IDs and the their time spent visiting the cue to the structure
            if strcmp(odorChecked, 'AIR')
                expData(fileIndex).listTotalIDsNearCueAIR=[tmpIDs1; tmpIDs2; tmpIDsInBoth]';
                expData(fileIndex).PrefIndexTotalIDsAIR=[tmpPIsOnly1; tmpPIsOnly2; tmpPIs]';
            elseif strcmp(odorChecked, 'CO2')           
                expData(fileIndex).listTotalIDsNearCueCO2=[tmpIDs1; tmpIDs2; tmpIDsInBoth]';
                expData(fileIndex).PrefIndexTotalIDsCO2=[tmpPIsOnly1; tmpPIsOnly2; tmpPIs]';
            elseif strcmp(odorChecked, 'postCO2')           
                expData(fileIndex).listTotalIDsNearCuePostCO2=[tmpIDs1; tmpIDs2; tmpIDsInBoth]';
                expData(fileIndex).PrefIndexTotalIDsPostCO2=[tmpPIsOnly1; tmpPIsOnly2; tmpPIs]';
            
            end
            
            % now Generate PI between base cue and random volume
            % check if any ID visited Base & random cues/vol
            if baseColorIndexList(fileIndex)== 1         
                %[tmpIDsInBoth, inTmpIDs1, inTmpIDs2]= intersect(tmpIDs1, tmpIDs2);
                [tmpIDsInBaseAndRnd, inTmpIDs1, inTmpIDs3]= intersect(tmpIDs1, tmpIDs3);
                %[tmpIDsInBoth, inBaseCueIDs, inTestCueIDs]= intersect(tmpIDs1, tmpIDs2);
            else
                %[tmpIDsInBoth, inTmpIDs1, inTmpIDs2]= intersect(tmpIDs1, tmpIDs2);
                [tmpIDsInBaseAndRnd, inTmpIDs2, inTmpIDs3]= intersect(tmpIDs2, tmpIDs3);
            end
            
            % Generate the PI values for Random Volume vs Base Cue           
            if any(tmpIDsInBaseAndRnd)
                % Assign PI values in function of where the Base cue was placed                
                if baseColorIndexList(fileIndex)== 1
                    %Base color was placed in position 1 (+Y axis) in current experiment
                    for i= 1:length(tmpIDsInBaseAndRnd)
                        tmpPIsInRnd(i,1)= (tmpCts3(inTmpIDs3(i))-tmpCts1(inTmpIDs1(i)))/(tmpCts3(inTmpIDs3(i))+tmpCts1(inTmpIDs1(i)));
                    end
                    % Remove ID visiting both cues as they are grouped in tmpIDsInBoth
                    % (Erase duplicate ID when it has visited both cues)
                    tmpIDsBase= setdiff(tmpIDs1, tmpIDsInBaseAndRnd);
                    tmpIDsRand= setdiff(tmpIDs3, tmpIDsInBaseAndRnd);
                else
                    %Base color was placed in position 2 (-Y axis) in current experiment
                    for i= 1:length(tmpIDsInBaseAndRnd)
                        tmpPIsInRnd(i,1)= (tmpCts3(inTmpIDs3(i))-tmpCts2(inTmpIDs2(i)))/(tmpCts3(inTmpIDs3(i))+tmpCts2(inTmpIDs2(i)));
                    end
                    % Remove ID visiting both cues as they are grouped in tmpIDsInBoth
                    % (Erase duplicate ID when it has visited both cues)
                    tmpIDsBase= setdiff(tmpIDs2, tmpIDsInBaseAndRnd);
                    tmpIDsRand= setdiff(tmpIDs3, tmpIDsInBaseAndRnd);
                 end
               
            else
                % ANY trajectory have visited BOTH cues.
                if baseColorIndexList(fileIndex)== 1
                    %Base color was placed in position 1 (+Y axis) in current experiment
                    tmpIDsBase= tmpIDs1;
                    tmpIDsRand= tmpIDs3;
                else
                    %Base color was placed in position 2 (-Y axis) in current experiment
                    tmpIDsBase= tmpIDs2;
                    tmpIDsRand= tmpIDs3;
                end
            end
            
            % Assign their Pref Index value for the trajectory that
            % only visited 1 of the volumes
            tmpPIsOnlyBase= ones(length(tmpIDsBase),1)*(-1);  %PI= -1
            tmpPIsOnlyRand= ones(length(tmpIDsRand),1);       %PI= 1

            %Add the IDs and the their time spent visiting the cue to the structure
            if strcmp(odorChecked, 'AIR')
                expData(fileIndex).listTotalIDsRndVolAIR=[tmpIDsBase; tmpIDsRand; tmpIDsInBaseAndRnd]';
                expData(fileIndex).PrefIndexTotalIDsRndVolAIR=[tmpPIsOnlyBase; tmpPIsOnlyRand; tmpPIsInRnd]';   
            elseif strcmp(odorChecked, 'CO2')           
                expData(fileIndex).listTotalIDsRndVolCO2=[tmpIDsBase; tmpIDsRand; tmpIDsInBaseAndRnd]';
                expData(fileIndex).PrefIndexTotalIDsRndVolCO2=[tmpPIsOnlyBase; tmpPIsOnlyRand; tmpPIsInRnd]'; 
            elseif strcmp(odorChecked, 'postCO2')           
                expData(fileIndex).listTotalIDsRndVolPostCO2=[tmpIDsBase; tmpIDsRand; tmpIDsInBaseAndRnd]';
                expData(fileIndex).PrefIndexTotalIDsRndVolPostCO2=[tmpPIsOnlyBase; tmpPIsOnlyRand; tmpPIsInRnd]';             
            end
        end
    end
end
clear tmpIDs1 tmpIDs2 tmpIDsBoth tmpPIsOnly1 tmpPIsOnly2 tmpPIs tmpCts1 tmpCts2
clear tmpIDsInBoth inTmpIDs1 inTmpIDs2  
clear sortedData split i fileIndex fileName filesList 
% ===========================
% ===========================


 

%% (13.) =================== Generate a SCATTER PLOT for color tested ===================
meanAllValues= zeros(length(testColorList),2);
meanValues= zeros(length(testColorList),1);
for typeIndex=1:length(testColorList)
    figure(gcf)
    % Find which entries in expData are related to each Tested color
    % expIndex give position regarding the cIndexes used to generate expData
    if workWithAllColors
        expIndex= find(strcmp(smryCounts.testedColor, testColorList(typeIndex)));
    else
        expIndex= find(strcmp(smryCounts.testedColor(cIndexes), testColorList(typeIndex)));
    end
    %PIsPerTypeCO2=[expData(expIndex).PrefIndexBothCuesCO2];
    allPIsPerTypeAIR= [expData(expIndex).PrefIndexTotalIDsAIR];
    allPIsPerTypeRndVolAIR= [expData(expIndex).PrefIndexTotalIDsRndVolAIR];
    allPIsPerTypeCO2= [expData(expIndex).PrefIndexTotalIDsCO2];
    allPIsPerTypeRndVolCO2= [expData(expIndex).PrefIndexTotalIDsRndVolCO2];
    %meanValues(typeIndex)= mean(PIsPerTypeCO2);
    meanAllValuesAIR(typeIndex,1)= mean(allPIsPerTypeAIR);
    meanAllValuesAIR(typeIndex,2)= mean(allPIsPerTypeRndVolAIR);
    meanAllValues(typeIndex,1)= mean(allPIsPerTypeCO2);
    meanAllValues(typeIndex,2)= mean(allPIsPerTypeRndVolCO2);
    %subplot(2,3,typeIndex)
    scatter(ones(size(allPIsPerTypeAIR,2),1), allPIsPerTypeAIR, 10, 'jitter', 'on', 'MarkerFaceColor', [0, 0.5, 0], 'MarkerFaceAlpha', 0.5);
    hold on
    scatter(ones(size(allPIsPerTypeRndVolAIR,2),1)*2, allPIsPerTypeRndVolAIR, 10, 'jitter', 'on', 'MarkerFaceColor', [0, 0, 0.5], 'MarkerFaceAlpha', 0.5);
    scatter(ones(size(allPIsPerTypeCO2,2),1)*3, allPIsPerTypeCO2, 10, 'jitter', 'on', 'MarkerFaceColor', [0, 0.5, 0], 'MarkerFaceAlpha', 0.5);
    scatter(ones(size(allPIsPerTypeRndVolCO2,2),1)*4, allPIsPerTypeRndVolCO2, 10, 'jitter', 'on', 'MarkerFaceColor', [0, 0, 0.5], 'MarkerFaceAlpha', 0.5);
    % plot also the means values 
    y= [meanAllValuesAIR(typeIndex,1), meanAllValuesAIR(typeIndex,1)];
    plot([0.9, 1.1], y, 'g', 'LineWidth',5,'MarkerEdgeColor',[0 0.5 0]);
    y= [meanAllValuesAIR(typeIndex,2), meanAllValuesAIR(typeIndex,2)];
    plot([1.9, 2.1], y, 'b', 'LineWidth',5,'MarkerEdgeColor',[0 0 0.5]);
    y= [meanAllValues(typeIndex,1), meanAllValues(typeIndex,1)];
    plot([2.9, 3.1], y, 'g', 'LineWidth',5,'MarkerEdgeColor',[0 0.5 0]);
    y= [meanAllValues(typeIndex,2), meanAllValues(typeIndex,2)];
    plot([3.9, 4.1], y, 'b', 'LineWidth',5,'MarkerEdgeColor',[0 0 0.5]);
    hold off
    %lg= [{' With AIR'}, {'With CO2'}];
    lg= [{'Base vs Test AIR'} {'Base vs Random AIR'} {'Base vs Test CO2'} {'Base vs Random CO2'}];
    legend(lg);
    title(strcat('PI near Tested color (',testColorList(typeIndex),') and h= 16cm (AIR and CO2'));
    xlim([0,5]);
    ylim([-1.2, 1.2]);
    ylabel('PI')
    if length(testColorList) > 1
        % Wait for user to plot data for new color
        disp('Pres "Enter to plot next Color');
        pause;
    end
    %outputFolder= 'AllColors_trajPI_scatterPlot\';
    %imgTitle= testColorList(typeIndex);
    %save_plot_in_exp_folder(gcf, imgTitle);
    figure(2);
    subplot(2,2,1);
    histogram(allPIsPerTypeAIR, 200);
    title(' Histogram PI values Test vs Base in AIR (h= 16cm)')
    subplot(2,2,2);
    histogram(allPIsPerTypeRndVolAIR, 20);
    title(' Histogram PI values Random vs Base in AIR (h= 16cm)')
    subplot(2,2,3);
    histogram(allPIsPerTypeCO2, 200);
    title(' Histogram PI values Test vs Base in CO2 (h= 16cm)')
    subplot(2,2,4);
    histogram(allPIsPerTypeRndVolCO2, 200);
    title(' Histogram PI values Random vs Base in CO2 (h= 16cm)')
end
clear typeIndex y imgTitle
% ===============================
% ===============================

%% (14.) =========== SAVE TRAJ. IDs PI VALUES PER TEST COLOR   =================
%   SAVE PI values per color used as tested cue
for c= unique(testColorList)'
    cIndexes= strcmp(cellstr(smryCounts.testedColor), char(c));
    testColorDataGrp= expData(cIndexes);
    %varName= strcat(char(c),'_PI_PerTraj_withRndVol_ForPaper.mat');
    varName= strcat(char(c),'_PI_PerTraj_With_h016.mat');
    if contains(outputFolder, 'AnStephensi')
        varPath= strcat(outputPath, 'AllColors_AnSteph_TrajPIs_withRndVol\');
    elseif contains(outputFolder, 'CxQuinquefasciatus')
        varPath= strcat(outputPath, 'AllColors_CxQuinq_TrajPIs_withRndVol\');
    else
        varPath= strcat(outputPath, 'AllColors_ForPaper_TrajPIs_withRndVol\');
    end
    save(strcat(varPath,varName), 'testColorDataGrp');
end
% ===============================
% ===============================


%% (15.) ============ VISITS IN VOLUMES OVER TIME===================
% - (WORKING WITH INDIVIDUAL .xls FILES PER EXPERIMENT)
% - Calculate # of trajectories visiting each volume over groups of time of X seconds
%   We will count how many visit occurred in each subset of time 
numGrps= 60;
odorChecked='postCO2'; %AIR - CO2 - postCO2
workWithAllColors= false;
% If workWithAllColors== false, choose the color to work with
colorGroup= 'GwT1';%'gray4.0'; %'R-HUE';   %Select the color to work with
%kk= filesList([1:5,10,14,16,20])

% Select with which data do you want to work (1 color/group or all colors) 
if workWithAllColors
    baseColorIndexList= smryCounts.baseColor_Pos;
    testColorList= unique(smryCounts.testedColor);
    % Load all files from group folder
    %filesPath= strcat(outputPath, 'AnalysisData_All_Exp_BaseC_Red\');
    %filesPath= strcat(outputPath, 'AnalysisData_WithRndVol_All_Exp_baseC_rhue\');
    filesPath= strcat(outputPath, 'AnalysisData_WithRndVol_All_Exp_baseC_white\');
    % Load the files in the individual folder
    cd(filesPath);
    filesList=dir(strcat('*countsInsideCueVol_',odorChecked,'.xlsx'));

else
    if contains(outputFolder, 'AnStephensi') || contains(outputFolder, 'CxQuinquefasciatus')
        baseColorIndexList= smryCounts.baseColor_Pos;
        testColorList= unique(smryCounts.testedColor);
        % Load the files in the individual folder
        filesPath= strcat(outputPath, outputFolder,'analysisData_newVersion\');
        % Load the files in the individual folder
        cd(filesPath);
        filesList=dir(strcat('*countsInsideCueVol_',odorChecked,'.xlsx'));
    else
        cIndexes= find(strncmp(smryCounts.testedColor, colorGroup, length(colorGroup)));
        baseColorIndexList= smryCounts.baseColor_Pos(cIndexes);
        testColorList= unique(smryCounts.testedColor(cIndexes));
        filesPath= strcat(outputPath, outputFolder,'analysisData_newVersion\');      
        %Build the names of the files to load
        filesList= cell2struct(strcat(smryCounts.date(cIndexes),'_countsInsideCueVol_',odorChecked,'.xlsx'), 'name', size(cIndexes,1));
        
        % Load  ALL the files in the individual folder
        %cd(filesPath);
        %filesList=dir(strcat('*countsInsideCueVol_',odorChecked,'.xlsx'));
    end
end

cd(workspace);
% generate how many visit in each volume per time group
[filesName, p1, p2, p3, p4]= load_ts_insect_in_volume_per_groups(numGrps, filesPath, filesList);
% group the visits in position per positions of the Base and Test cues
for fileIndex= 1:length(filesName)
    if workWithAllColors
        % Find the corresponding index in the summary report
        smryIndex= find(strcmp(smryCounts.date, filesName(1)));
    else
        smryIndex= fileIndex;
    end
    
    if baseColorIndexList(smryIndex) == 1
        visitBaseC= p1(fileIndex,:);
        visitTestC= p2(fileIndex,:);
    else
        visitBaseC= p2(fileIndex,:);
        visitTestC= p1(fileIndex,:);
    end
    
   
    % Group the information in the visitOverTime structure
    switch(odorChecked)
        case 'AIR'
            % Set file Name
            visitOverTime(fileIndex).date= filesName(fileIndex);
            visitOverTime(fileIndex).nearTestAIR=visitTestC;
            visitOverTime(fileIndex).nearBaseAIR=visitBaseC;
            visitOverTime(fileIndex).nearRandAIR= p3(fileIndex,:);
        case 'CO2'          
            visitOverTime(fileIndex).nearTestCO2=visitTestC;
            visitOverTime(fileIndex).nearBaseCO2=visitBaseC;
            visitOverTime(fileIndex).nearRandCO2= p3(fileIndex,:);
        case 'postCO2'          
            visitOverTime(fileIndex).nearTestPostCO2=visitTestC;
            visitOverTime(fileIndex).nearBasePostCO2=visitBaseC;
            visitOverTime(fileIndex).nearRandPostCO2= p3(fileIndex,:);
    end
    clear visitBaseC visitTestC smryIndex
end
clear fileIndex
% ===================================================
% ==================================


%% (16.)========== PLOT VISIT NEAR CUES OVER TIME =========
% (Need to the improved/finished to plot for all colors)
colors={'red', 'blue', 'black', 'green', 'yellow', 'red', 'blue', 'black', 'green', 'yellow', 'red', 'blue'};
figure()

if workWithAllColors
    % Generate plots for each of the color used
    for colorGroup=testColorList
        cIndexes= find(strncmp(smryCounts.testedColor, colorGroup, 4));
        % Generate mean values for the group of experiment sharing same test cue 
        subplot(2,2,1)
        hold on
        plot(smooth([visitOverTime(i).nearTestAIR, visitOverTime(i).nearTestCO2]), colors{i})
        plot([numGrps, numGrps], [0, 100], ':')
        hold off
        subplot(2,2,2)
        hold on
        plot(cumsum([visitOverTime(i).nearTestAIR, visitOverTime(i).nearTestCO2]), colors{i})
        plot([numGrps, numGrps], [0, 100], ':')
        hold off
        subplot(2,2,3)
        hold on
        plot(smooth([visitOverTime(i).nearBaseAIR, visitOverTime(i).nearBaseCO2]), colors{i})
        plot([numGrps, numGrps], [0, 100], ':')
        hold off
        subplot(2,2,4)
        hold on
        plot(smooth([visitOverTime(i).nearRandAIR, visitOverTime(i).nearRandCO2]), colors{i})
        plot([numGrps, numGrps], [0, 100], ':')
        hold off    
    end
    subplot(2,2,1)
    title('Visit over time Test cue AIR-CO2')
    subplot(2,2,2)
    title('cumsum Visit over time Test cue AIR-CO2')
    subplot(2,2,3)
    title('Visit over time Base cue AIR-CO2')
    subplot(2,2,4)
    title('Visit over time Random Vol AIR-CO2')
else
    % Generate the plots for the given color only
    % Generate mean values for the group of experiment sharing same test cue 
    subplot(2,2,1)
    hold on
    plot(smooth([visitOverTime(i).nearTestAIR, visitOverTime(i).nearTestCO2]), testColorList)
    plot([numGrps, numGrps], [0, 100], ':')
    hold off
    subplot(2,2,2)
    hold on
    plot(cumsum([visitOverTime(i).nearTestAIR, visitOverTime(i).nearTestCO2]), testColorList)
    plot([numGrps, numGrps], [0, 100], ':')
    hold off
    subplot(2,2,3)
    hold on
    plot(smooth([visitOverTime(i).nearBaseAIR, visitOverTime(i).nearBaseCO2]), testColorList)
    plot([numGrps, numGrps], [0, 100], ':')
    hold off
    subplot(2,2,4)
    hold on
    plot(smooth([visitOverTime(i).nearRandAIR, visitOverTime(i).nearRandCO2]), testColorList)
    plot([numGrps, numGrps], [0, 100], ':')
    hold off    
    subplot(2,2,1)
    title('Visit over time Test cue AIR-CO2')
    subplot(2,2,2)
    title('cumsum Visit over time Test cue AIR-CO2')
    subplot(2,2,3)
    title('Visit over time Base cue AIR-CO2')
    subplot(2,2,4)
    title('Visit over time Random Vol AIR-CO2')
end

for i=1:12
    subplot(2,2,1)
    hold on
    plot(smooth([visitOverTime(i).nearTestAIR, visitOverTime(i).nearTestCO2]), colors{i})
    plot([numGrps, numGrps], [0, 100], ':')
    hold off
    subplot(2,2,2)
    hold on
    plot(cumsum([visitOverTime(i).nearTestAIR, visitOverTime(i).nearTestCO2]), colors{i})
    plot([numGrps, numGrps], [0, 100], ':')
    hold off
    subplot(2,2,3)
    hold on
    plot(smooth([visitOverTime(i).nearBaseAIR, visitOverTime(i).nearBaseCO2]), colors{i})
    plot([numGrps, numGrps], [0, 100], ':')
    hold off
    subplot(2,2,4)
    hold on
    plot(smooth([visitOverTime(i).nearRandAIR, visitOverTime(i).nearRandCO2]), colors{i})
    plot([numGrps, numGrps], [0, 100], ':')
    hold off    
end
subplot(2,2,1)
title('Visit over time Test cue AIR-CO2')
subplot(2,2,2)
title('cumsum Visit over time Test cue AIR-CO2')
subplot(2,2,3)
title('Visit over time Base cue AIR-CO2')
subplot(2,2,4)
title('Visit over time Random Vol AIR-CO2')


figure()
plot(smooth([visitOverTime(5).nearTestAIR, visitOverTime(5).nearTestCO2],0.25), 'red')
hold on
plot(smooth([visitOverTime(5).nearTestAIR, visitOverTime(5).nearTestCO2],0.95), 'black')
plot(smooth([visitOverTime(5).nearTestAIR, visitOverTime(5).nearTestCO2],0.1,'rloess'), 'blue')

hold off
% ===========================
% ===========================

%% (17.) ============ SAVE VISITS OVER TIME PER GROUP OF TIME (X sec duration) ================
%   SAVE visits near volumes over time values per color used as tested cue
%  
for c= unique(testColorList)'
    if exist('visitOverTime', 'var')
        if workWithAllColors
            cIndexes= strcmp(cellstr(smryCounts.testedColor), char(c));
            % We are working with visit over time groups analysis
            testColorDataGrp= visitOverTime(cIndexes);
        else
            testColorDataGrp= visitOverTime;
        end
    end
    % We are working with the overall PI per trajectories
%     if workWithAllColors
%         testColorDataGrp= expData(cIndexes);
%     else
%         testColorDataGrp= expData;
%     end
    varName= strcat(char(c),'_visitsNearCues_OverTime_withRndVol.mat');
    varPath= strcat(outputPath, 'visitNearCues_overTime_ForPaper_withRndVol\');
    save(strcat(varPath,varName), 'testColorDataGrp');
end
% ===============================
% ===============================



%% WARNING: In following analysis the postCO2 is not analyzed as we are not
% sure how much the CO2 remains can affect insect behavior


% (18) ====================== (II) % of trajectories NEAR TESTED CUE (1 color)
% If working the RAW data : TrajSummary contains information about the trajectories. Column values are: 
% [Number of Traj AIR, Avg duration traj AIR, Max. duration traj AIR, ...
%   Number of Traj CO2, Avg duration traj CO2, Max. duration traj CO2, ...   
%       Number of Traj PostCO2, Avg duration traj PostCO2, Max. duration traj PostCO2]
% If working with SUMMARY REPORT data, use the smryTrajectorys table

figure()
clear temp tempNorm clear trajNearTestedCue
tempAIR=[];
tempCO2=[];
for v=i'
    % group the data in function of the base color cue position
    if (smryCounts.baseColor_Pos(v) == 1)
        tempAIR= [tempAIR; p2IDsAIR(v)];
        tempCO2= [tempCO2; p2IDsCO2(v)];
%            tempPostCO2= [tempPostCO2; p2IDsPostCO2(v)];
    else
        tempAIR= [tempAIR; p1IDsAIR(v)];
        tempCO2= [tempCO2; p1IDsCO2(v)];
%            tempPostCO2= [tempPostCO2; p1IDsPostCO2(v)];
    end
end
% Calculate the % of trajectories interested in the tested color cue
tempNorm= bsxfun(@rdivide,(tempAIR*100),smryTrajectories.numOfTrajAIR(i));
trajNearTestedCue(:,1)= tempNorm;
%tempValListAIR(i)= tempNorm;
tempNorm= bsxfun(@rdivide,(tempCO2*100),smryTrajectories.numOfTrajCO2(i));
trajNearTestedCue(:, 2)= tempNorm;
%tempValListCO2(i)= tempNorm;
%    tempNorm= bsxfun(@rdivide,(tempPostCO2*100),smryTrajectories.numOfTrajPostCO2(i));
%    trajNearCues(:, 3)= tempNorm;
%    tempValListPostCO2(i)= tempNorm;

%Generate a SCATTER PLOT for this type of mutation
subplot(3,3,typeIndex)
scatter(ones(size(trajNearTestedCue,1),1), trajNearTestedCue(:,1), 'jitter', 'on', 'MarkerFaceColor', [0, 0.5, 0]);
hold on
scatter(ones(size(trajNearTestedCue,1),1), trajNearTestedCue(:,2), 'jitter', 'on', 'MarkerFaceColor', [0.5, 0, 0]);
%    scatter(ones(size(trajNearCues,1),1), trajNearCues(:,3), 'jitter', 'on', 'MarkerFaceColor', [0, 0, 0.5]);

% plot also the means values 
y= [mean(trajNearTestedCue(:,1)), mean(trajNearTestedCue(:,1))];
plot([0.9, 1.1], y, '--g', 'MarkerSize',15,'MarkerEdgeColor',[0 0.5 0]);
y= [mean(trajNearTestedCue(:,2)), mean(trajNearTestedCue(:,2))];
plot([0.9, 1.1], y, '--r', 'MarkerSize',15,'MarkerEdgeColor',[0.5 0 0]);
%    y= [mean(trajNearCues(:,3)), mean(trajNearCues(:,3))];
%    plot([0.9, 1.1], y, '--b', 'MarkerSize',15,'MarkerEdgeColor',[0 0 0.5]);
hold off
lg= [{' With AIR'}, {'With CO2'}];
legend(lg);
title(strcat('% of trajectories detected near black cue (',typeValue,')'));
xlim([0,2]);
ylim([0, 25]);
ylabel('% of Trajectories')
typeIndex= typeIndex+1;    

clear tempAIR tempCO2 tempPostCO2 tempNorm
clear lg y i typeValue tipeIndex
% ==================================================
% ==================================================    









%% == == Christian DIANGCO scripts == == 

% ============================
% Get trajectories that start before CO2 is active and end when CO2
% is active.
trajInBoth = cell(length(filesList),1);
for fileIndex= 1:length(filesList)
    fileName= filesList(fileIndex).name;
    disp(strcat(' - Working with file: ', {' '}, fileName));
    trajInBoth{fileIndex} = traj_in_prevCO2_and_withCO2(dataset(fileIndex),flightTimeLimit);
end

%Once you have the trajectories, let's do the following 
% - save the list of IDs for the trajectories found (one ID per trajectory)
%   and the date of the experiment (it is in the dataset.fileName) in an excel
%   file. Information from all experiments in the current dataset must be
%   saved in the same xlsx file and the Excel file should be saved in path:
%       C:\Users\Riffell Lab\Desktop\Christian_Diangco\matlab\output\exp_4_grays
% - Estimate the direction in the XY and XZ planes at any given time Ti. We
%   can consider it as vectors from position T(i-1) to position T(i) (for all
%   XYZ positions for this ID)
%       - put in excel file

% Create array of IDs in trajInBoth with their corresponding experiment
% dates and save to ids_in_prevCO2_and_withCO2.xlsx
ids_and_dates= zeros(1,2);
for fileIndex= 1:length(filesList)
    id = -1;
    currentData= trajInBoth{fileIndex,1};
    if ~isempty(currentData)
        uniqueIds= unique(cell2mat(currentData(:,1)));
        expDate= str2num(cell2mat(extractBetween(dataset(fileIndex).fileName,'','_')));
        for id= transpose(uniqueIds)
            ids_and_dates= cat(1,ids_and_dates,[id,expDate]);
        end
        continue;
    end
    % if there are no trajectories in both pre and during CO2 for the
    % experiment, the id in the excel file for that experiment will be -1
    expDate= str2num(cell2mat(extractBetween(dataset(fileIndex).fileName,'','_')));
    ids_and_dates= cat(1,ids_and_dates,[id,expDate]);
end
ids_and_dates(1,:)= [];
%filePath= 'C:\Users\Riffell Lab\Desktop\Christian_Diangco\matlab\output\exp_4_grays\ids_in_prevCO2_and_withCO2.xlsx';
filePath= 'C:\Users\chris\OneDrive\Documents\MATLAB\output\ids_in_prevCO2_and_withCO2.xlsx';
writematrix(ids_and_dates,filePath);

% Find estimates of direction in XY and XZ planes at each timestamp for
% each insect. Direction is the angle(degrees) of the vector from timestamp T-1 to
% T.

for fileIndex= 1:length(filesList)
    % looks at trajectories that start in prevCO2 and end in withCO2, but
    % should work for large dataset
    currentData= trajInBoth{fileIndex,1};
    if ~isempty(currentData)
        currentData= estimate_direction(currentData);

    % creates excel files for each date, with dxy and dxz appended to the
    % right side of the original array for that date
    %filePath= 'C:\Users\Riffell Lab\Desktop\Christian_Diangco\matlab\output\exp_4_grays\';
    filePath= 'C:\Users\chris\OneDrive\Documents\MATLAB\output\';
    expDate= cell2mat(extractBetween(dataset(fileIndex).fileName,'','_'));
    fileName= strcat(expDate,'_direction_at_ts_TEST','.xlsx');
    filePath = strcat(filePath,fileName);
    writecell(currentData, filePath);
    end
end

% == == END == ==




%% ======================================================================
% ==================================== Old Code ==============================
% ======================================================================
