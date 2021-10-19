

%function datasetList= load_all_dataset()
close all;
clear variables;
format short;
%opengl('save', 'software');

fileID= '20191230_110918';
fileName= strcat(fileID,'.mainbrain.h5');
%Define the base color to use in the analysis
baseColor= 'white'; % Must be black || white || R-HUE
plotDataFlag= false;

% Information related to the data to analyze
workspace= 'C:\Users\dalon\Documents\GitHub\AnalysisFLYDRA_2020';
inputPath= 'D:\DataFlydra\MosquitoProject\INPUT_Data\';
outputPath= 'D:\DataFlydra\MosquitoProject\OUTPUT_Graphs\';
expDep= 'expDeprecated\';

subFolder= 'FLYDRA_Trajectories_mutants_lineX\';

%Create the outputFolder folder
outputFolder = load_output_folder(outputPath, subFolder);

% Specify the file to work with 
inputPath= strcat(inputPath,subFolder);


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


% Load file name to work with 
filePath= strcat(inputPath,fileName);
% Load all the information from the FLYDRA .h5 file    
[attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z]= load_data_from_file(filePath, loadFullDataset);
%tempData= [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z];

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
indexWithCO2= find(dataset.attr_time(:) >= ts_startCO2 & dataset.attr_time(:) < ts_endCO2-1);
indexPostCO2= find(dataset.attr_time(:) >= ts_endCO2);

% Define which odor stime was active at each timestamp in teh dataset
dataset.stim(indexPrevCO2,1)= {'AIR'};
dataset.stim(indexWithCO2,1)= {'CO2'};
dataset.stim(indexPostCO2,1)= {'postCO2'};

%Load the position where the baseCue is and the color tested in experiment
baseColorIndex=find(strcmp(cuesSetup([2,3],1), baseColor));
testedColor= cuesSetup(find(~strcmp(cuesSetup([2,3],1), baseColor))+1,1);

%Clear temporary variables from workspace
clear fileName loadFullDataset 
clear attr_id attr_time attr_frame attr_x attr_y attr_z
clear  cuesSetup ts_startAIR ts_startCO2 ts_endCO2 ts_endAIR mType mGender



% ===========================    
% Count the amount of trajectories detected while stim is prior, during and
% post CO2
trajSummary= zeros(length(dataset), 9);
%Find the indexes fro AIR, CO2 and postCO2

% Count the total number of trajectories longer than flightTimeLimit
% seconds. AS each row value in dataset.attr_id(:) is related to 1
% frame we only need to check in how many frames that ID has been
% detected

for i=1:3
    if i == 1
        data= [dataset.attr_id(indexPrevCO2) dataset.attr_time(indexPrevCO2)];
        disp(' --- Data prev CO2 --');  
    elseif i==2 
        data= [dataset.attr_id(indexWithCO2) dataset.attr_time(indexWithCO2)];
        disp(' --- Data with CO2 --');  
    else
        data= [dataset.attr_id(indexPostCO2) dataset.attr_time(indexPostCO2)];
        disp(' --- Data post CO2 --');  

    end
    data= [dataset.attr_id() dataset.attr_time()];
    uniqueID=unique(data(:,1));
    %Transpose from a column matrix to a row matrix
    uniqueID=uniqueID';
    % Initialice local variables
    totalTraj=0; 
    flightTimeList= nan(length(uniqueID),1);
    totalFlightTime=0;
    % For each insect ID, check if 
    for index= 1:length(uniqueID)
        objID= uniqueID(index);
        %Load the frames where appears the current objID
        objTime= find(data(:,1) == objID);
        % Estimate the duration of each insect ID trajectory
        duration= get_trajectory_duration(data(objTime(:),2));

        %Consider only if flight duration is bigger than flightTimeLimit seconds
        if duration >= flightTimeLimit
            totalTraj= totalTraj+1;
            %Keep temporary stored all trajectories duration
            flightTimeList(index)=duration;
            totalFlightTime= totalFlightTime+duration;
        end
    end
    % Calculate average and max trajectory duration
    avgTime= totalFlightTime/totalTraj;
    maxTime= max(flightTimeList);
    disp(strcat(' * Total amount of trajectories over ',num2str(flightTimeLimit),' seconds: ',num2str(totalTraj)));
    %disp(strcat(' * Average flight time for this dataset (seconds): ',num2str(totalFlightTime/length(uniqueID))));
    disp(strcat(' * Average flight time for this dataset (seconds): ',num2str(avgTime)));
    disp(strcat(' * Shortest flight duration (seconds): ',num2str(min(flightTimeList))));
    disp(strcat(' * Longest flight duration (seconds): ',num2str(max(maxTime))));
    disp(' ----- ');
        if i == 1
        data= [dataset.attr_id(indexPrevCO2) dataset.attr_time(indexPrevCO2)];
        disp(' --- Data prev CO2 --');  
    elseif i==2 
        data= [dataset.attr_id(indexWithCO2) dataset.attr_time(indexWithCO2)];
        disp(' --- Data with CO2 --');  
    else
        data= [dataset.attr_id(indexPostCO2) dataset.attr_time(indexPostCO2)];
        disp(' --- Data post CO2 --');  

    end 
    
end
disp(' --- Data prev CO2 --');  
[trajSummary(1,1), trajSummary(1,2), trajSummary(1,3)]= count_trajectories([dataset.attr_id(indexPrevCO2) dataset.attr_time(indexPrevCO2)], flightTimeLimit);
disp(' --- Data during CO2 --');
[trajSummary(1,4), trajSummary(1,5), trajSummary(1,6)]= count_trajectories([dataset.attr_id(indexWithCO2), dataset.attr_time(indexWithCO2)], flightTimeLimit);
disp(' --- Data post CO2 --');
[trajSummary(1,7), trajSummary(1,8), trajSummary(1,9)]= count_trajectories([dataset.attr_id(indexPostCO2), dataset.attr_time(indexPostCO2)], flightTimeLimit);
disp(' =====================================');








%+++++++++++++++++++++++
kkPath= 'D:\DataFlydra\MosquitoProject\OUTPUT_Graphs\exp_mutants_lineX\analysisData_oldVersion\';
fileName= '20191226_113013_countsInsideCueVol_AIR.xlsx';
disp(strcat(' - Working with file: ', {' '}, fileName));
% data loaded structure: [vCuePosition	objID	timeStamp NaN initialTimeStamp finalTimeStamp]
dataFromExcel = table2array(readtable(strcat(kkPath, fileName)));

size(unique(dataFromExcel(:,2)));

disp(strcat(' - Working with file: ', {' '}, fileName));
    filesName(fileIndex)= {fileName(1:8)};
    % data loaded structure: [vCuePosition	objID	timeStamp NaN initialTimeStamp finalTimeStamp]
    dataFromExcel = table2array(readtable(strcat(filesPath, fileName)));
    expData(fileIndex).name= fileName;
    % Initialize the correct fields in function of the odor used
    if strcmp(odorChecked, 'AIR')
        expData(fileIndex).listIDsInCueAIR=[];
        expData(fileIndex).listTimeInCueAIR=[];
        expData(fileIndex).listTrajTimeAIR=[];
        expDataAllTrajSz(fileIndex).listIDsInCueAIR=[];
        expDataAllTrajSz(fileIndex).listTimeInCueAIR=[];
        expDataAllTrajSz(fileIndex).listTrajTimeAIR=[];
    elseif strcmp(odorChecked, 'CO2')
        expData(fileIndex).listIDsInCueCO2=[];
        expData(fileIndex).listTimeInCueCO2=[];
        expData(fileIndex).listTrajTimeCO2=[];
        expDataAllTrajSz(fileIndex).listIDsInCueCO2=[];
        expDataAllTrajSz(fileIndex).listTimeInCueCO2=[];
        expDataAllTrajSz(fileIndex).listTrajTimeCO2=[];
    end
    if isempty(dataFromExcel)
        disp(strcat('   -> No counts detected in file: ',{' '}, fileName));
    else
        %Initialize the object that will contain the information for the current experiment            
        tmpListTime=[];
        tmpIndex= 1;
        % To load timestamp data (time spent by insect near thew cue
        % Sorted the the data by cue visited ID timestamp
        sortedData= sortrows(dataFromExcel(:,1:3));
        % Find the first appearence of the position 2 
        split= find(sortedData(:,1)==2, 1);
        if baseColorIndexList(expIndex)== 1
            %Base color was placed in position 1 (+Y axis) in current experiment
            sortedData= sortedData(split:end,:);
        else
            %Base color was placed in position 2 (-Y axis) in current experiment
            sortedData= sortedData(1:(split-1),:);
        end       
        %Estimates the time spend for each ID in the pos where was TESTED cue
        tmpIndex=1;
        for objID= unique(sortedData(:,2))'
            idIndexes= find(sortedData(:,2) ==objID);
            % Calculate the Delta time between timestamps values for the id 
            tsDiff= diff(sortedData(idIndexes,3));
            % find which  sequential counts didn't happen "~=consecutively"
            out= find(tsDiff > 0.05);
            t=0;
            k=1;
            if any(out)
                for i= out
                    % the insect left the volume at a given moment
                    t= t + sum(tsDiff(k:(i-1)));
                    k= i+1;
                end
                t= t + sum(tsDiff(k:end));
            else
                %the insects has been inside the volume all the time
                t= sum(tsDiff);
            end
            tmpListTime(1,tmpIndex)= t;
            tmpIndex= tmpIndex+1;
        end
        %Add the IDs and the their time spent visiting the cue to the structure
        if strcmp(odorChecked, 'AIR')
            expData(fileIndex).listIDsInCueAIR=unique(sortedData(:,2))';
            expData(fileIndex).listTimeInCueAIR=tmpListTime;
            expDataAllTrajSz(fileIndex).listIDsInCueAIR=unique(sortedData(:,2))';
            expDataAllTrajSz(fileIndex).listTimeInCueAIR=tmpListTime;     
        elseif strcmp(odorChecked, 'CO2')
            expData(fileIndex).listIDsInCueCO2=unique(sortedData(:,2))';
            expData(fileIndex).listTimeInCueCO2=tmpListTime;
            expDataAllTrajSz(fileIndex).listIDsInCueCO2=unique(sortedData(:,2))';
            expDataAllTrajSz(fileIndex).listTimeInCueCO2=tmpListTime;
        end
    end









