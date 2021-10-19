% Function to load the timestamp values and order and positions of clues
% (odor/color and their XYZ position) inside the test section for the mosquito project
% - Arguments:
%       - fileName: Name of the current file that need its settings to be uploaded in the variables workspace
% - Return (a row matrix with):
%       - cuesSetup: XYZ positions for the odor and the visual cues specified as:
%           {'Odor': xOdor yOdor zOdorZ; 'color1' x1 y1 z1; 'color2' x2 y2 z2; ... 'colorN' xN yN zN}
%       - ts_startAIR: timestamp showing the begining of the experiment
%       - ts_startCO2: timestamp showing the begining of the CO2 release
%       - ts_endCO2: timestamp showing the end of the CO2 release
%       - ts_endCO2: timestamp showing the end of the experiment
%       - mosquitoType: type of mosquitoes used in the experiment: wt (wild type), m1 (mutant 1), m2 (mutant 2), m0 (mutant 1 and 2), 
%           lx (x= 1..N) values are assigned in the experiments that we run blinded (UCSB knows which type are but us not)
%       - mosquitoGender: gender of the mosquitoes used in the experiment: female (f) and male(m)
function [cuesSetup, ts_startAIR, ts_startCO2, ts_endCO2, ts_endAIR, mosquitoType, mosquitoGender]=  load_exp_settings(fileName)

    switch(fileName)
          
% * WHITE AND RED PERPENDICULAR TO AIRFLOW --------------------------                 
        case '20200629_081404.mainbrain.h5'
            mosquitoType= 'wt';
            mosquitoGender= 'f';           
            cuesSetup={'odor' 0.15 0 0.20; 'white' 0.48 0.09 0.01; 'RS1' 0.48 -0.09 0.01};
            % timeStamp value when the MasFlowController with AIR was initialize
            ts_startAIR= 1593443664; 
            % timeStamp value when the MasFlowController with CO2 was initialize
            ts_startCO2= 1593447264;
            % timeStamp value when the MasFlowController with CO2 was stopped
            ts_endCO2= 1593450864;
            % timestamp value when the MassFlowController wtih CO2 was stopped
            % (END OF EXPERIMENT)
            ts_endAIR= 1593454465;   
     
% * WHITE AND GRAY PERPENDICULAR TO AIRFLOW --------------------------            
        case '20200218_104800.mainbrain.h5'
            mosquitoType= 'wt';
            mosquitoGender= 'f';           
            cuesSetup={'odor' 0.15 0 0.20; 'gray9.5' 0.47 0.07 0.01; 'white' 0.47 -0.07 0.01};
            % timeStamp value when the MasFlowController with AIR was initialize
            ts_startAIR= 1582051770; 
            % timeStamp value when the MasFlowController with CO2 was initialize
            ts_startCO2= 1582055370;
            % timeStamp value when the MasFlowController with CO2 was stopped
            ts_endCO2= 1582062570;
            % timestamp value when the MassFlowController wtih CO2 was stopped
            % (END OF EXPERIMENT)
            ts_endAIR= 1582069770;
       
            
        otherwise
            disp(strcat('  * ERROR: Settings for experiment',{' '}, fileName, {' '},' not added yet'))
            disp('        Please add the settings for this experiment in file load_exp_settings.m') 
            disp('  * Check if settings are set in script: load_exp_settings_all_exp_types.m');
            mosquitoType= 'None';
            mosquitoGender= 'None';
            cuesSetup={'odor' 0.15 0 0.20; 'black' 0.995 -0.14 -0.2; 'white' 0.995 0.12 -0.2};
            %cuesSetup={'odor' 0.15 0 0.20; 'black' 0.32 0.14 0.02; 'white' 0.31 -0.12 0.02};
            
            % timeStamp value when the MasFlowController with AIR was initialize
            ts_startAIR=0;
            % timeStamp value when the MasFlowController with CO2 was initialize
            ts_startCO2=0;
            % timeStamp value when the MasFlowController with CO2 was stopped
            ts_endCO2=0;
            % timestamp value when the MassFlowController wtih CO2 was stopped
            % (END OF EXPERIMENT)
            ts_endAIR=0;
    end

end