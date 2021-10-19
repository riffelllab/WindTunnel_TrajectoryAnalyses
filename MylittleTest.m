clear all;
format long;

index= 1;
fps=60.0;
% Lower time threshold (in seconds) to plot a trajectory or not
flightTimeLimit= 1.5;

% Test section X,Y,Z limits (start_x= -lim_x and end_x= +lim_x --> Total x distance is lim_x*2)
% These lim values are the current size of the Test Section (DO NOT CHANGE IT)
lim_x= 0.9144;
lim_y= 0.3048;
lim_z= 0.6096;

% Load file name to work with 
fileName= 'MyLittleTest';

% Load all the information from the FLYDRA .h5 file    
%[attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z]= load_data_from_file(filePath, loadFullDataset);
sz= 100;

attr_id= sort(randi(sz/5, sz,1));
attr_time= [1578683325.44615;1578683325.45438;1578683325.46956;1578683325.48428;1578683325.49502;1578683325.51159;1578683325.53435;1578683325.54475;1578683325.56591;1578683325.57967;1578683325.59680;1578683325.60911;1578683325.62633;1578683325.64392;1578683325.66125;1578683325.67683;1578683325.69411;1578683325.71084;1578683325.72820;1578683325.74501;1578683325.76204;1578683325.78171;1578683325.79681;1578683325.81029;1578683325.82810;1578683325.85147;1578683325.86157;1578683325.88210;1578683325.89365;1578683325.91612;1578683325.92610;1578683325.94975;1578683325.96566;1578683325.98277;1578683325.99929;1578683326.01075;1578683326.03309;1578683326.04398;1578683326.06205;1578683326.07872;1578683326.16989;0;1578683326.18454;1578683326.19483;1578683326.20863;1578683326.22324;1578683326.23727;1578683326.25002;1578683326.27818;1578683326.29302;1578683326.30532;1578683326.31716;1578683326.33068;1578683326.18229;1578683326.19271;1578683326.20624;1578683326.22087;1578683326.23444;1578683326.24778;1578683326.27370;1578683326.28926;1578683326.30277;1578683326.31485;1578683326.32877;1578683325.43120;1578683325.43758;1578683325.45171;1578683325.46784;1578683325.48208;1578683325.49239;1578683325.50916;1578683325.53214;1578683325.54248;1578683325.56422;1578683325.57808;1578683325.59379;1578683325.60747;1578683325.62444;1578683325.64173;1578683325.65834;1578683325.67482;1578683325.69169;1578683325.70825;1578683325.72544;1578683325.74217;1578683325.75931;1578683325.77823;1578683325.79342;1578683325.80782;1578683325.82515;1578683325.84619;1578683325.85851;1578683325.87742;1578683325.89082;1578683325.91345;1578683325.92372;1578683325.94707;1578683325.96286;1578683325.97879;1578683325.99665;1578683326.00825];
%attr_frame=[219;220;221;222;223;224;225;226;227;228;229;230;231;232;233;234;235;236;237;238;239;240;241;242;243;244;245;246;247;248;249;250;251;252;253;254;255;256;257;258;261;262;263;264;265;266;267;268;269;270;271;272;273;263;264;265;266;267;268;269;270;271;272;273;218;219;220;221;222;223;224;225;226;227;228;229;230;231;232;233;234;235;236;237;238;239;240;241;242;243;244;245;246;247;248;249;250;251;252;253]
attr_frame= sort(190+randperm(100)');
attr_x= (lim_x -(-lim_x)).*rand(sz,1) + (-lim_x);
attr_y= (lim_y -(-lim_y)).*rand(sz,1) + (-lim_y);
attr_z= (lim_z - 0).*rand(sz,1) + 0;

% Load the different stimuli used and their position in the test
% section
%[cuesSetup,ts_startAIR,ts_startCO2,ts_endCO2,ts_endAIR, mType, mGender]=  load_exp_settings(fileName);
mType= 'wt';
mGender= 'f';
cuesSetup={'odor' 0.15 0 0.20; 'red' 0.45 0.10 0.01; 'white' 0.65 0.10 0.01; 'black' 0.45 -0.10 0.01; 'green' 0.65 -0.10 0.01};
% timeStamp value when the MasFlowController with AIR was initialize
ts_startAIR=attr_time(1);
% timeStamp value when the MasFlowController with CO2 was initialize
ts_startCO2=attr_time(round(length(attr_time)/3));
% timeStamp value when the MasFlowController with CO2 was stopped
ts_endCO2=attr_time(round(length(attr_time)/3 + length(attr_time)/3));
% timestamp value when the MassFlowController wtih CO2 was stopped
% (END OF EXPERIMENT)
ts_endAIR=attr_time(end);
% Empty the entries recorded before and after the mass flow
% controller script is launched.
expIndexes= find(attr_time(:) >= ts_startAIR & attr_time(:) <= ts_endAIR);

% Estimate the timestamps values for odor stimulus ON and OFF
% attr_time contains timestamp epochs (in seconds)
indexPrevCO2= find(attr_time(:) < ts_startCO2);
indexWithCO2= find(attr_time >= ts_startCO2 & attr_time(:) < ts_endCO2-1);
indexPostCO2= find(attr_time >= ts_endCO2);

%Fill the dataset with the data from current experiment
dataset(index).fileName=    fileName;
dataset(index).type=        mType;
dataset(index).gender=       mGender;
dataset(index).expCues=    cuesSetup;
dataset(index).attr_id=     attr_id;
dataset(index).attr_time=   attr_time;
dataset(index).attr_frame=  attr_frame;
dataset(index).attr_x=      attr_x;
dataset(index).attr_y=      attr_y;
dataset(index).attr_z=      attr_z;
dataset(index).stim(indexPrevCO2,1)= {'AIR'};
dataset(index).stim(indexWithCO2,1)= {'CO2'};
dataset(index).stim(indexPostCO2,1)= {'postCO2'};

plot_XY_XZ_heatmaps_v2(dataset(index).attr_x, dataset(index).attr_y, dataset(index).attr_z, dataset(index).expCues, 'myLittleTest')