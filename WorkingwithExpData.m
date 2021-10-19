clear all; close all; 
dataWhite = load('C:\Users\Claire Rusch\Downloads\expData_mutants_nearWhiteCue.mat');
dataBlack = load('C:\Users\Claire Rusch\Downloads\expData_mutants_dataset.mat');

%% All trajectories > 1.5 sec
cd('C:/Users/Claire Rusch/Desktop')

for i = 1:length(dataWhite.expData)
    % Percent Time In Cue
    PercentAirWhite = dataWhite.expData(i).percentNearCueAIR;
    PercentAirBlack = dataBlack.expData(i).percentNearCueAIR;
    PercentCWhite = dataWhite.expData(i).percentNearCueCO2;
    PercentCBlack = dataBlack.expData(i).percentNearCueCO2;
    % Total Time Trajectory    
    TimeAirWhite = dataWhite.expData(i).listTrajTimeAIR;
    TimeAirBlack = dataBlack.expData(i).listTrajTimeAIR;
    TimeCWhite = dataWhite.expData(i).listTrajTimeCO2;
    TimeCBlack = dataBlack.expData(i).listTrajTimeCO2;
    % Create Table    
    PercentAirBlack = table(PercentAirBlack(PercentAirBlack~= 0)', TimeAirBlack(PercentAirBlack ~= 0)', 'VariableNames', {'PercentCue', 'TotalTime'});
    PercentAirWhite = table(PercentAirWhite(PercentAirWhite~=0)', TimeAirWhite(PercentAirWhite ~= 0)', 'VariableNames', {'PercentCue', 'TotalTime'});
    PercentCBlack = table(PercentCBlack(PercentCBlack ~= 0)', TimeCBlack(PercentCBlack ~= 0)', 'VariableNames', {'PercentCue', 'TotalTime'});
    PercentCWhite = table(PercentCWhite(PercentCWhite ~= 0)', TimeCWhite(PercentCWhite ~= 0)', 'VariableNames', {'PercentCue', 'TotalTime'});
    
    %writetable(PercentAirWhite,'PercentWhitecue_Air.xls','Sheet', dataWhite.expData(i).name(1:8))
    %writetable(PercentCWhite,'PercentWhitecue_CO2.xls','Sheet',dataWhite.expData(i).name(1:8))
    %writetable(PercentAirBlack,'PercentBlackcue_Air.xls','Sheet',dataBlack.expData(i).name(1:8))
    %writetable(PercentCBlack,'PercentBlackcue_CO2.xls','Sheet',dataBlack.expData(i).name(1:8))
end

%% Only trajectories that connected 1 cue. 

for i = 1:length(dataWhite.expData)
    % Find common ID
    IDBlackAir = dataBlack.expData(i).listIDsInCueAIR;
    IDWhiteAir = dataWhite.expData(i).listIDsInCueAIR;
    IDBlackCO2 = dataBlack.expData(i).listIDsInCueCO2;
    IDWhiteCO2 = dataWhite.expData(i).listIDsInCueCO2;
    CommonIDAIR = intersect(dataBlack.expData(i).listIDsInCueAIR, dataWhite.expData(i).listIDsInCueAIR);
    CommonIDCO2 = intersect(dataBlack.expData(i).listIDsInCueCO2, dataWhite.expData(i).listIDsInCueCO2);
    % Indices common ID
    IDBlackAir = find(~ismember(IDBlackAir,CommonIDAIR));
    IDWhiteAir = find(~ismember(IDWhiteAir,CommonIDAIR));
    IDBlackCO2 = find(~ismember(IDBlackCO2,CommonIDCO2));
    IDWhiteCO2 = find(~ismember(IDWhiteCO2,CommonIDCO2));
    % Percent Time In Cue
    PercentAirWhite = dataWhite.expData(i).percentNearCueAIR(IDWhiteAir);
    PercentAirBlack = dataBlack.expData(i).percentNearCueAIR(IDBlackAir);
    PercentCWhite = dataWhite.expData(i).percentNearCueCO2(IDWhiteCO2);
    PercentCBlack = dataBlack.expData(i).percentNearCueCO2(IDBlackCO2);
    % Total Time Trajectory    
    TimeAirWhite = dataWhite.expData(i).listTrajTimeAIR(IDWhiteAir);
    TimeAirBlack = dataBlack.expData(i).listTrajTimeAIR(IDBlackAir);
    TimeCWhite = dataWhite.expData(i).listTrajTimeCO2(IDWhiteCO2);
    TimeCBlack = dataBlack.expData(i).listTrajTimeCO2(IDBlackCO2);
    % Create Table    
    PercentAirBlack = table(PercentAirBlack(PercentAirBlack~= 0)', TimeAirBlack(PercentAirBlack ~= 0)', 'VariableNames', {'PercentCue', 'TotalTime'});
    PercentAirWhite = table(PercentAirWhite(PercentAirWhite~=0)', TimeAirWhite(PercentAirWhite ~= 0)', 'VariableNames', {'PercentCue', 'TotalTime'});
    PercentCBlack = table(PercentCBlack(PercentCBlack ~= 0)', TimeCBlack(PercentCBlack ~= 0)', 'VariableNames', {'PercentCue', 'TotalTime'});
    PercentCWhite = table(PercentCWhite(PercentCWhite ~= 0)', TimeCWhite(PercentCWhite ~= 0)', 'VariableNames', {'PercentCue', 'TotalTime'});
    
    %writetable(PercentAirWhite,'PercentWhitecue_Air.xls','Sheet', dataWhite.expData(i).name(1:8))
    %writetable(PercentCWhite,'PercentWhitecue_CO2.xls','Sheet',dataWhite.expData(i).name(1:8))
    %writetable(PercentAirBlack,'PercentBlackcue_Air.xls','Sheet',dataBlack.expData(i).name(1:8))
    %writetable(PercentCBlack,'PercentBlackcue_CO2.xls','Sheet',dataBlack.expData(i).name(1:8))
    
    % Time on cue Trajectory    
    TimeAirWhite = dataWhite.expData(i).listTimeInCueAIR(IDWhiteAir);
    TimeAirBlack = dataBlack.expData(i).listTimeInCueAIR(IDBlackAir);
    TimeCWhite = dataWhite.expData(i).listTimeInCueCO2(IDWhiteCO2);
    TimeCBlack = dataBlack.expData(i).listTimeInCueCO2(IDBlackCO2);
    
    
    PICountsAir(i) = (length(PercentAirBlack.PercentCue) - length(PercentAirWhite.PercentCue))/(length(PercentAirBlack.PercentCue) + length(PercentAirWhite.PercentCue));
    PITimeAir(i) = (sum(TimeAirBlack) - sum(TimeAirWhite))/ (sum(TimeAirBlack) + sum(TimeAirWhite));
    PICountsCO2(i) = (length(PercentCBlack.PercentCue) - length(PercentCWhite.PercentCue))/(length(PercentCBlack.PercentCue) + length(PercentCWhite.PercentCue));
    PITimeCO2(i) = (sum(TimeCBlack) - sum(TimeCWhite))/ (sum(TimeCBlack) + sum(TimeCWhite));
    date{i} = dataWhite.expData(i).name(1:8);
    %writetable(PercentAirWhite,'PercentWhitecue_Air.xls','Sheet', dataWhite.expData(i).name(1:8))
    %writetable(PercentCWhite,'PercentWhitecue_CO2.xls','Sheet',dataWhite.expData(i).name(1:8))
    %writetable(PercentAirBlack,'PercentBlackcue_Air.xls','Sheet',dataBlack.expData(i).name(1:8))
    %writetable(PercentCBlack,'PercentBlackcue_CO2.xls','Sheet',dataBlack.expData(i).name(1:8))
end

T = table(date', PICountsAir', PITimeAir', PICountsCO2', PITimeCO2', 'VariableNames', {'Date', 'PIcountsAir', 'PItimeAir', 'PIcountsCO2', 'PItimeCO2'});
writetable(T, 'PImutants.xls', 'Sheet', 1)


%% Only trajectory that connected bobth cue to compute a PI (sigh)
%% Only trajectories that connected 1 cue. 

for i = 1:length(dataWhite.expData)
    % Find common ID
    IDBlackAir = dataBlack.expData(i).listIDsInCueAIR;
    IDWhiteAir = dataWhite.expData(i).listIDsInCueAIR;
    IDBlackCO2 = dataBlack.expData(i).listIDsInCueCO2;
    IDWhiteCO2 = dataWhite.expData(i).listIDsInCueCO2;
    CommonIDAIR = intersect(dataBlack.expData(i).listIDsInCueAIR, dataWhite.expData(i).listIDsInCueAIR);
    CommonIDCO2 = intersect(dataBlack.expData(i).listIDsInCueCO2, dataWhite.expData(i).listIDsInCueCO2);
    % Indices common ID
    IDBlackAir = find(ismember(IDBlackAir,CommonIDAIR));
    IDWhiteAir = find(ismember(IDWhiteAir,CommonIDAIR));
    IDBlackCO2 = find(ismember(IDBlackCO2,CommonIDCO2))
    IDWhiteCO2 = find(ismember(IDWhiteCO2,CommonIDCO2))
    
    % Percent Time In Cue
    PercentAirWhite = dataWhite.expData(i).percentNearCueAIR(IDWhiteAir);
    PercentAirBlack = dataBlack.expData(i).percentNearCueAIR(IDBlackAir);
    PercentCWhite = dataWhite.expData(i).percentNearCueCO2(IDWhiteCO2);
    PercentCBlack = dataBlack.expData(i).percentNearCueCO2(IDBlackCO2);
    % Total Time Trajectory    
    TimeAirWhite = dataWhite.expData(i).listTrajTimeAIR(IDWhiteAir);
    TimeAirBlack = dataBlack.expData(i).listTrajTimeAIR(IDBlackAir);
    TimeCWhite = dataWhite.expData(i).listTrajTimeCO2(IDWhiteCO2);
    TimeCBlack = dataBlack.expData(i).listTrajTimeCO2(IDBlackCO2);
    % Create Table    
    PercentAirBlack = table(PercentAirBlack(PercentAirBlack~= 0)', TimeAirBlack(PercentAirBlack ~= 0)', 'VariableNames', {'PercentCue', 'TotalTime'});
    PercentAirWhite = table(PercentAirWhite(PercentAirWhite~=0)', TimeAirWhite(PercentAirWhite ~= 0)', 'VariableNames', {'PercentCue', 'TotalTime'});
    PercentCBlack = table(PercentCBlack(PercentCBlack ~= 0)', TimeCBlack(PercentCBlack ~= 0)', 'VariableNames', {'PercentCue', 'TotalTime'});
    PercentCWhite = table(PercentCWhite(PercentCWhite ~= 0)', TimeCWhite(PercentCWhite ~= 0)', 'VariableNames', {'PercentCue', 'TotalTime'});
    
    %writetable(PercentAirWhite,'PercentWhitecue_Air.xls','Sheet', dataWhite.expData(i).name(1:8))
    %writetable(PercentCWhite,'PercentWhitecue_CO2.xls','Sheet',dataWhite.expData(i).name(1:8))
    %writetable(PercentAirBlack,'PercentBlackcue_Air.xls','Sheet',dataBlack.expData(i).name(1:8))
    %writetable(PercentCBlack,'PercentBlackcue_CO2.xls','Sheet',dataBlack.expData(i).name(1:8))
    
    % Time on cue Trajectory    
    TimeAirWhite = dataWhite.expData(i).listTimeInCueAIR(IDWhiteAir);
    TimeAirBlack = dataBlack.expData(i).listTimeInCueAIR(IDBlackAir);
    TimeCWhite = dataWhite.expData(i).listTimeInCueCO2(IDWhiteCO2);
    TimeCBlack = dataBlack.expData(i).listTimeInCueCO2(IDBlackCO2);
    PIAir = [];
    PICO2 = [];
    for j = 1:length(TimeAirWhite)
        PIAir(j) = (TimeAirBlack(j)-TimeAirWhite(j))/(TimeAirBlack(j)+TimeAirWhite(j));
    end
    for j = 1:length(TimeCWhite)
        PICO2(j) = (TimeCBlack(j) - TimeCWhite(j))/(TimeCBlack(j)+TimeCWhite(j));
    end
    if length(TimeAirWhite) > 0
        T = table(PIAir', 'VariableNames', {'PI'});
        writetable(T,'IndividualPI_Air.xls','Sheet', dataWhite.expData(i).name(1:8))
    end
    if length(TimeCWhite)>0 
        T2 = table(PICO2', 'VariableNames', {'PI'});
        writetable(T2,'IndividualPI_CO2.xls','Sheet',dataWhite.expData(i).name(1:8))
    end
    %PICountsAir(i) = (length(PercentAirBlack.PercentCue) - length(PercentAirWhite.PercentCue))/(length(PercentAirBlack.PercentCue) + length(PercentAirWhite.PercentCue));
    %PITimeAir(i) = (sum(TimeAirBlack) - sum(TimeAirWhite))/ (sum(TimeAirBlack) + sum(TimeAirWhite));
    %PICountsCO2(i) = (length(PercentCBlack.PercentCue) - length(PercentCWhite.PercentCue))/(length(PercentCBlack.PercentCue) + length(PercentCWhite.PercentCue));
    %PITimeCO2(i) = (sum(TimeCBlack) - sum(TimeCWhite))/ (sum(TimeCBlack) + sum(TimeCWhite));
    %date{i} = dataWhite.expData(i).name(1:8);
    %writetable(PercentAirBlack,'PercentBlackcue_Air.xls','Sheet',dataBlack.expData(i).name(1:8))
    %writetable(PercentCBlack,'PercentBlackcue_CO2.xls','Sheet',dataBlack.expData(i).name(1:8))
end
