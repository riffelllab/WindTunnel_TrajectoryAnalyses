%Function to plpot preference plots for each of the axis of a volume (X,Y,Z
%axis)
function check_preference_in_testSection(data, imgTitle, fullExpDatasets)

    
    if fullExpDatasets
        posX= data(:,1);
        posY= data(:,2);
        posZ= data(:,3);
    else
        posX= data(:,4);
        posY= data(:,5);
        posZ= data(:,6);
    end;
    
    meanX= mean(posX);
    meanY= mean(posY);
    meanZ= mean(posZ);
    
    switch (imgTitle(end-6:end))
    case 'prevCO2'
        t2= '(Intial air only)';
    case 'withCO2'
        t2= '(Air + CO2)';
    case 'postCO2'
        t2= '(Final air only)';
    otherwise
        if imgTitle(3)== '_'
            % type and gender of mosquito provided in imgTitle
            t2= strcat(imgTitle(1:2),'-', imgTitle(4), '-','dataset');
        else
            t2= '(Full dataset)';
        end;
    end;
    figTitle= sprintf(strcat('Distribution of each axis in TS','\n',t2));
    
    
    figure()
    subplot(1,3,1);
    % X axis histogram
    hist(posX);
    hold on;
    line([meanX, meanX], ylim, 'LineWidth', 2, 'Color', 'r', 'LineStyle',':');    
    hold off;
    xlabel('Test Section X axis');
    ylabel('Counts')
    subplot(1,3,2);
    % Y axis histogram
    hist(posY);
    hold on;
    line([meanY, meanY], ylim, 'LineWidth', 2, 'Color', 'r', 'LineStyle',':');    
    hold off;
    xlabel('Test Section Y axis');
    ylabel('Counts')
    title(figTitle); 
    subplot(1,3,3);
    % Z axis histogram
    hist(posZ);
    hold on;
    line([meanZ, meanZ], ylim, 'LineWidth', 2, 'Color', 'r', 'LineStyle',':');
    hold off;
    xlabel('Test Section Z axis');
    ylabel('Counts')

    print(gcf,strcat(imgTitle, '.png'),'-dpng','-r300');        % *// 300 dpi
end