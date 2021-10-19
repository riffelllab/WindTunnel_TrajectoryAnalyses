%function that plot the main information from the Summary report
% Arguments:
%   - smryTrajectories: table with the Number of trajectories counted,
%       their avg time duration & the longest trajectorie detected grouped in:
%       prevCO2, WithCO2 and PostCO2)
%   - smryCounts table with the counts per volume information for each experiment
%               stored in the Summary report
%   - smryPIs: table with the preference index values for each experiment
%               stored in the Summary report
%   - expList: List with the type of experiments loaded from the summary
%               report
function plot_summaries(smryTrajectories, smryCounts, smryPIs, expList)

    %Load the output path
    outputPath = evalin('base', 'outputPath');

    
    % ==== Plot Pref Indexes ====
    if ~empty(smryPIs)
        data= str2double(string(table2array(smryPIs(:,2))));
        bar(data);
        %xAxis= 'Experiments',
        %set(gca,'XtickLabel', xAxis);
        if contains(version, 'R2019')
            xtickangle(45)    
        end
        title('Preference Index towards the tested color');
        print(gcf,strcat(outputPath,  't_prefIndexSmryPlot', '.png'),'-dpng','-r300');        % *// 300 dpi
      
    
        % ==== Plot Preference indexes grouped per color & with their +/- std error ====
        numOfExp= [6, 5, 3, 3, 3, 3];
        data= str2double(string(table2array(smryPIs(:,2))));
        blackRed= data(1:6);
        whiteRed= data(7:11);
        whiteGray1= [data(12), data(16), data(23)];
        whiteGray2= [data(17),data(18), data(20)];
        whiteGray3= [data(13), data(16), data(21)];
        whiteGray4= [data(14), data(19), data(22)];
    
        % Estimates they mean value for each color 
        meanPIperColor= [mean(blackRed), mean(whiteRed), mean(whiteGray1), mean(whiteGray2), mean(whiteGray3), mean(whiteGray4)];
        xAxis= { 'Black vs Red', 'White vs Red', 'White vs Gray1', 'White vs Gray2', 'White vs Gray3', 'White vs Gray4'};
        % Estimate the error for the mean values
        errList= [std(blackRed)/numOfExp(1), std(whiteRed)/numOfExp(2), std(whiteGray1)/numOfExp(3), std(whiteGray2)/numOfExp(4), std(whiteGray1)/numOfExp(5), std(whiteGray1)/numOfExp(6)];
    
        % Sort the Data & Rearrange Labels
        [sortedPIs, newIndices] = sort(meanPIperColor); % sorts in *ascending* order
        sortedLabels = xAxis(newIndices); 
        sortedErr= errList(newIndices);
    
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
    
        set(gca,'XtickLabel', sortedLabels, 'xtick',1:length(sortedPIs));
    
        if contains(version, 'R2019')
            xtickangle(45)    
        end
        title('Preference Index towards the tested color');
        print(gcf,strcat(outputPath,  't_PIgroupedSmryPlot', '.png'),'-dpng','-r300');        % *// 300 dpi
    end        
    
        
    % ==== Plot the counts inside each volume surrounding the visual cue ====
    if ~empty(smryCounts)
        data= str2double(string(table2array(smryCounts(:,2:3))));
        bar(data);
        %xAxis= 'Experiments',
        %set(gca,'XtickLabel', xAxis);
        if contains(version, 'R2019')
            xtickangle(45)    
        end
        title('Counts per Position');
        print(gcf,strcat(outputPath,  't_countInVolSmryPlot', '.png'),'-dpng','-r300');        % *// 300 dpi  
    end    
    

    % ==== Plot the information about the trajectories ====
    if ~empty(smryTrajectories)
        % Load the data in 2 separate arrays
        prevCO2= str2double(string(table2array(smryTrajectories(:,2:4))));
        withCO2= str2double(string(table2array(smryTrajectories(:,5:7))));
        postCO2= str2double(string(table2array(smryTrajectories(:,8:10))));

        %MEAN values trajectories
        trjMeanList= [mean(prevCO2(1:6,1)) mean(prevCO2(7:11,1)) mean(prevCO2(12:end,1)); ...
                      mean(withCO2(1:6,1)) mean(withCO2(7:11,1)) mean(withCO2(12:end,1)); ...
                      mean(postCO2(1:6,1)) mean(postCO2(7:11,1)) mean(postCO2(12:end,1))];

        %MEAN values average time of flight
        avgTimeMeanList= [mean(prevCO2(1:6,2)) mean(prevCO2(7:11,2)) mean(prevCO2(12:end,2)); ...
                           mean(withCO2(1:6,2)) mean(withCO2(7:11,2)) mean(withCO2(12:end,2)); ...
                           mean(postCO2(1:6,2)) mean(postCO2(7:11,2)) mean(postCO2(12:end,2))];          

        %MEAN values trajectories
        maxTimeMeanList= [mean(prevCO2(1:6,3)) mean(prevCO2(7:11,3)) mean(prevCO2(12:end,2)); ...
                           mean(withCO2(1:6,3)) mean(withCO2(7:11,3)) mean(withCO2(12:end,2)); ...
                           mean(postCO2(1:6,3)) mean(postCO2(7:11,3)) mean(postCO2(12:end,2))];          
    
        subplot(2,1,1);
        bar(trjMeanList);
        xAxis= {'PrevCO2', 'WithCO2', 'PostCO2'};
        set(gca,'XtickLabel', xAxis);
        if contains(version, 'R2019')
            xtickangle(45)    
        end
        title('mean trajectories counted (per Exp.)');
        legend(expList);
    

        subplot(2,1,2);
        bar(avgTimeMeanList);
        % hold on;
        % bar(maxTimeMeanList);
        % hold off;
        xAxis= {'PrevCO2', 'WithCO2', 'PostCO2'};
        set(gca,'XtickLabel', xAxis);
        if contains(version, 'R2019')
            xtickangle(45)    
        end
        title('average trajectories duration (per Exp.)');
        legend(expList);
    
        print(gcf,strcat(outputPath,  't_trajSmryPlot', '.png'),'-dpng','-r300');        % *// 300 dpi
    end
end