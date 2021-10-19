% ===========================
% Load instants where insect are inside a cue volume and group them by time after the odor was released
odorChecked='CO2';
if (exist('outputPath', 'var') &&  exist('outputFolder', 'var'))
    % Load files
    filesPath= strcat(outputPath, outputFolder);
    cd(filesPath);
    filesList=dir(strcat('*countsInsideCueVol_',odorChecked,'.xlsx'));
    cd(workspace);
    % Define number of groups to generate
    numGrps=120; %=120 =60 sec per group || =4 =1800 sec per grp
    %find the counts near visual cue over time and per experiment
    %[filesName, p1, p2, p3, p4]= load_ts_insect_in_volume_per_groups(numGrps, filesPath, filesList);
    [filesName, p1, p2, p3, p4, t1, t2]= load_insect_data_per_time_groups(numGrps, filesPath, filesList);
    % save this information in a .CSV file in the project folder
    save_counts_over_time_in_file(p1, p2, p3, p4)
else
    disp('  * Warning: file path and/or folder not specified in the variabes workspace ');
end
%If loading AIR or postCO2 data, assign it to its correct variable
switch(odorChecked)
    case'CO2'
        p1CO2= p1;
        p2CO2= p2;
        filesName= filesName';
        t1CO2= t1;
        t2CO2= t2;
    case'AIR'
        p1AIR= p1;
        p2AIR= p2;
        filesNameAIR= [filesName(1,1:5),{'20200225'},{'20200226'},filesName(1,6:end)]';
        t1AIR= t1;
        t2AIR=t2;
        %Add missing files to dataset
        p1AIR= [p1AIR(1:5,:); zeros(2,length(p1AIR(1,:))); p1AIR(6:end,:)];
        p2AIR= [p2AIR(1:5,:); zeros(2,length(p2AIR(1,:))); p2AIR(6:end,:)];
        t1AIR= [t1AIR(1:5,:); zeros(2,length(p1AIR(1,:))); t1AIR(6:end,:)];
        t2AIR= [t2AIR(1:5,:); zeros(2,length(p2AIR(1,:))); t2AIR(6:end,:)];
    case'postCO2'
        p1PostCO2= p1;
        p2PostCO2= p2;
        filesNamePostCO2= filesName';
        t1PostCO2= t1;
        t2PostCO2= t2;
end
%If you are only using 2 Visual Cues
clear p1 p2 p3 p4 t1 t2


% =============================
%As experiment after 2020.03.23 is 3 hours duration (instead of 5), withCO2 and
%postCO2 are only 1 hours long. group the time groups per couples (t and
%t+1) to make each time group == 1 minute -> filesName(13)== '20200331'
for expIndex=13:length(filesName)
    p1CO2(expIndex,:)= horzcat(p1CO2(expIndex,1:2:end)+p1CO2(expIndex,2:2:end), zeros(1,60));
    p2CO2(expIndex,:)= horzcat(p2CO2(expIndex,1:2:end)+p2CO2(expIndex,2:2:end), zeros(1,60));
    t1CO2(expIndex,:)= horzcat(t1CO2(expIndex,1:2:end)+t1CO2(expIndex,2:2:end), zeros(1,60));
    t2CO2(expIndex,:)= horzcat(t2CO2(expIndex,1:2:end)+t2CO2(expIndex,2:2:end), zeros(1,60));

    p1PostCO2(expIndex,:)= horzcat(p1PostCO2(expIndex,1:2:end)+p1PostCO2(expIndex,2:2:end), zeros(1,60));
    p2PostCO2(expIndex,:)= horzcat(p2PostCO2(expIndex,1:2:end)+p2PostCO2(expIndex,2:2:end), zeros(1,60));
    t1PostCO2(expIndex,:)= horzcat(t1PostCO2(expIndex,1:2:end)+t1PostCO2(expIndex,2:2:end), zeros(1,60));
    t2PostCO2(expIndex,:)= horzcat(t2PostCO2(expIndex,1:2:end)+t2PostCO2(expIndex,2:2:end), zeros(1,60));


end
% The aclimatation part is ALWAYS 1 hour duration (for all exp types)
p1AIR= horzcat(p1AIR(:,1:2:end)+p1AIR(:,2:2:end));
p2AIR= horzcat(p2AIR(:,1:2:end)+p2AIR(:,2:2:end));
t1AIR= horzcat(t1AIR(:,1:2:end)+t1AIR(:,2:2:end));
t2AIR= horzcat(t2AIR(:,1:2:end)+t2AIR(:,2:2:end));




% ================ BOXPLOT FOR GRAY VS WHITE =================
% * Do the number of trajectories recorded on cues change between experiments?
indexG1= find(strcmp(testedColorList(:),'gray1'));
indexG2= find(strcmp(testedColorList(:),'gray2'));
indexG3= find(strcmp(testedColorList(:),'gray3'));
indexG4= find(strcmp(testedColorList(:),'gray4'));

boxplot([trajSummary(indexG1,4), trajSummary(indexG2,4),trajSummary(indexG3,4), trajSummary(indexG4,4)])
title('Number of trajectories detected when CO2 is ON per visual cue')
xNames=[{'Gray1'}, {'Gray2'},{'Gray3'}, {'Gray4'}];
set(gca,'XtickLabel', xNames);
ylabel('Trajectories detected')
% ==================================================
% ==================================================    


% ==================================================
% * Do the dynamics change between experiments
% ==================================================
%load the type of grays
colorCues= unique(testedColorList)';
% Initialize variable to obtain the mean and std error per colorCue  
i=1;
duration=60;
meanList= zeros(length(colorCues), duration);
errList= zeros(length(colorCues), duration);
meanListPostCO2= zeros(length(colorCues), duration);
errListPostCO2= zeros(length(colorCues), duration);
meanListAIR= zeros(length(colorCues), duration);
errListAIR= zeros(length(colorCues), duration);
meanTimeListAIR= zeros(length(colorCues), duration);
errTimeListAIR= zeros(length(colorCues), duration);
meanTimeListCO2= zeros(length(colorCues), duration);
errTimeListCO2= zeros(length(colorCues), duration);
meanTimeListPostCO2= zeros(length(colorCues), duration);
errTimeListPostCO2= zeros(length(colorCues), duration);
for c=colorCues
    indexes= find(strcmp(testedColorList, c));
    %Sum the counts per position for each of the time groups
    sumPerColor=cumsum(p1CO2(indexes,1:duration)+p2CO2(indexes,1:duration), 2);
    sumPerColorPostCO2=cumsum(p1PostCO2(indexes,1:duration)+p2PostCO2(indexes,1:duration),2);
    sumPerColorAIR= cumsum(p1AIR(indexes,1:duration)+p2AIR(indexes,1:duration),2);
    % Sum the times spend by some insects in the volume, is there any
    % difference with the counts plot?
    sumTimeInPosCO2= cumsum(t1CO2(indexes, 1:duration)+t2CO2(indexes, 1:duration), 2);
    sumTimeInPosAIR= cumsum(t1AIR(indexes, 1:duration)+t2AIR(indexes, 1:duration), 2);
    sumTimeInPosPostCO2= cumsum(t1PostCO2(indexes, 1:duration)+t2PostCO2(indexes, 1:duration), 2);
    if length(indexes) > 1          
        %If several exp, generates the mean over their cummulative sum over time
        tempNorm= bsxfun(@rdivide,(sumPerColor*100),sumPerColor(:,end));
        tempNorm(isnan(tempNorm))=0;
        meanPerColor= mean(tempNorm);
        errCO2= std(tempNorm)/sqrt(length(indexes));

        tempNorm= bsxfun(@rdivide,(sumPerColorAIR*100),sumPerColor(:,end));
        tempNorm(isnan(tempNorm))=0;
        meanPerColorAIR= mean(tempNorm);
        errAIR= std(tempNorm)/sqrt(length(indexes));

        tempNorm= bsxfun(@rdivide,(sumPerColorPostCO2*100),sumPerColor(:,end));
        tempNorm(isnan(tempNorm))=0;
        meanPerColorPostCO2= mean(tempNorm);
        errPostCO2= std(tempNorm)/sqrt(length(indexes));

        %Working with the time matrices
        tempNorm= bsxfun(@rdivide,(sumTimeInPosAIR*100),sumTimeInPosCO2(:,end));
        tempNorm(isnan(tempNorm))=0;
        meanTimeAIR= mean(tempNorm);
        errTimeAIR= std(tempNorm)/sqrt(length(indexes));
        tempNorm= bsxfun(@rdivide,(sumTimeInPosCO2*100),sumTimeInPosCO2(:,end));
        tempNorm(isnan(tempNorm))=0;
        meanTimeCO2= mean(tempNorm);
        errTimeCO2= std(tempNorm)/sqrt(length(indexes));
        tempNorm= bsxfun(@rdivide,(sumTimeInPosPostCO2*100),sumTimeInPosCO2(:,end));
        tempNorm(isnan(tempNorm))=0;
        meanTimePostCO2= mean(tempNorm);
        errTimePostCO2= std(tempNorm)/sqrt(length(indexes));
    end

    meanList(i,:)= meanPerColor;
    errList(i,:)= errCO2;
    meanListPostCO2(i,:)= meanPerColorPostCO2;
    errListPostCO2(i,:)= errPostCO2;
    meanListAIR(i,:)= meanPerColorAIR;
    errListAIR(i,:)= errAIR;
    meanTimeListCO2(i,:)= meanTimeCO2;
    errTimeListCO2(i,:)= errTimeCO2;
    meanTimeListAIR(i,:)= meanTimeAIR;
    errTimeListAIR(i,:)= errTimeAIR;
    meanTimeListPostCO2(i,:)= meanTimePostCO2;
    errTimeListPostCO2(i,:)= errTimePostCO2;
    i=i+1;
end
% SHADED ERROR BAR
x= 1:duration; %[1:60; 1:60; 1:60; 1:60]
%Plot the counts withCO2 vs postCO2 each plot is for a value of gray
figure(1);
figure(2);
title('Counts in volume withCO2 vs postCO2 (per color)');
%set(gcf,'Title', title);
lg= [{'Pre CO2'}, {'With CO2'},{'Post CO2'}];
for i= 1:4
    figure(1);
    pp1=subplot(2,2,i);
    shadedErrorBar(x, meanListAIR(i,:)', errListAIR(i,:)','lineprops','-b')
    %pp1.YLim= [0 max(max(meanList))+2000];
    pp1.YLim= [0 150];
    hold on;
    shadedErrorBar(x, meanList(i,:)', errList(i,:)','lineprops','-r')
    shadedErrorBar(x, meanListPostCO2(i,:)', errListPostCO2(i,:)','lineprops','-g')
    hold off;
    xlabel('Minute #');
    %xlim([1, 120]);
    ylabel(' % Counted ');
    title(colorCues(i));
    legend(lg, 'Location','eastoutside');

    % Second figure showing time 
    figure(2)
    pp2=subplot(2,2,i);
    shadedErrorBar(x, meanTimeListAIR(i,:)', errTimeListAIR(i,:)','lineprops','-b')
    %pp1.YLim= [0 max(max(meanList))+2000];
    pp2.YLim= [0 150];
    hold on;
    shadedErrorBar(x, meanTimeListCO2(i,:)', errTimeListCO2(i,:)','lineprops','-r')
    shadedErrorBar(x, meanTimeListPostCO2(i,:)', errTimeListPostCO2(i,:)','lineprops','-g')
    hold off;
    xlabel('Minute #');
    %xlim([1, 120]);
    ylabel(' % of time  ');
    title(colorCues(i));
    legend(lg, 'Location','eastoutside');

end   
% ==================================================
% ==================================================    


% ==================================================
% Plot the cumSum of the 3 parts of the experiment per 1 type of gray
% ==================================================
indexes= find(strcmp(testedColorList, 'gray4'));
sumPerColor=cumsum(p1CO2(indexes,1:duration)+p2CO2(indexes,1:duration), 2);
sumPerColorPostCO2=cumsum(p1PostCO2(indexes,1:duration)+p2PostCO2(indexes,1:duration),2);
sumPerColorAIR= cumsum(p1AIR(indexes,1:duration)+p2AIR(indexes,1:duration),2);
maxCts= max(max(max(sumPerColorAIR(:,end)), max(sumPerColor(:,end))), max(sumPerColorPostCO2(:,end)));

% This plot requires plots extra infromation as data from the
% relFlightActDuration matrix
figure()
title('Counts in volume withCO2 vs postCO2 (per color)');
%set(gcf,'Title', title);
lg= [{'Pre CO2'}, {'With CO2'},{'Post CO2'}];
for i= 1:length(indexes)
    pp1=subplot(2,3,i);
    plot(sumPerColorAIR(i, :)');
    hold on
    plot(sumPerColor(i,:)');
    plot(sumPerColorPostCO2(i,:)');
    hold off
    xlabel('Minute #');
    ylim([1, maxCts]);
    ylabel(' Counts ');
    title(filesName(indexes(i)));
    legend(lg, 'Location','eastoutside');
end
pp1= subplot(2,3, 6);
bp= bar(relFlightActDuration(indexes,:));
l= cell(1,3);
l{1}='Prior CO2'; l{2}='With CO2'; l{3}= 'Post CO2';
legend(bp,l, 'Location', 'northeast');
set(gca,'XtickLabel', xNames(indexes));
if contains(version, 'R2019')
    xtickangle(90)    
end
%title('Estimated flight activity compared to prior-CO2 activity');
title('Estimated flight activity compared to previous activity');   
% ==================================================
% ==================================================


% ==================================================
% Analyze behavior for a type of GRAY
% ==================================================
% 0. See plot estimatedActivity_vs_prevCO2 or
% estActivity_vs_prevActivity to detect weird experiments
% 1. plot all the AIR/CO2/postCO2 shaded error plots (mean and std)
% 2. plot  data without the days that behave weird behavior
% 3. plot the experiments with weird behavior (more activity in AIR or/and postCO2 than with CO2  

x= 1:duration;
grayX= 4; %index for GRAY3 in the mean and strd error lists

% Exp indexes good and weird 
weirdIndex= [8];
goodIndex= [3, 11, 17, 18];

figure()
%Plot All exp for GRAYX
pp1=subplot(2,2,1);
shadedErrorBar(x, meanListAIR(grayX,:)', errListAIR(grayX,:)','lineprops','-b')
%pp1.YLim= [0 max(max(meanList))+2000];
pp1.YLim= [0 150];
hold on;
shadedErrorBar(x, meanList(grayX,:)', errList(grayX,:)','lineprops','-r')
shadedErrorBar(x, meanListPostCO2(grayX,:)', errListPostCO2(grayX,:)','lineprops','-g')
hold off;
xlabel('Minute #');
%xlim([1, 120]);
ylabel(' % Counted ');
title(colorCues(grayX));
legend(lg, 'Location','eastoutside');

% 2. plot  data without the days that behave weird behavior  
sumPerColor=cumsum(p1(goodIndex,1:duration)+p2(goodIndex,1:duration), 2);
sumPerColorPostCO2=cumsum(p1PostCO2(goodIndex,1:duration)+p2PostCO2(goodIndex,1:duration),2);
sumPerColorAIR= cumsum(p1AIR(goodIndex,1:duration)+p2AIR(goodIndex,1:duration),2);
if length(goodIndex) > 1          
    %If several exp, generates the mean over their cummulative sum over time
    tempNorm= bsxfun(@rdivide,(sumPerColor*100),sumPerColor(:,end));
    tempNorm(isnan(tempNorm))=0;
    meanPerColor= mean(tempNorm);
    errCO2= std(tempNorm)/sqrt(length(goodIndex));

    tempNorm= bsxfun(@rdivide,(sumPerColorAIR*100),sumPerColor(:,end));
    tempNorm(isnan(tempNorm))=0;
    meanPerColorAIR= mean(tempNorm);
    errAIR= std(tempNorm)/sqrt(length(goodIndex));

    tempNorm= bsxfun(@rdivide,(sumPerColorPostCO2*100),sumPerColor(:,end));
    tempNorm(isnan(tempNorm))=0;
    meanPerColorPostCO2= mean(tempNorm);
    errPostCO2= std(tempNorm)/sqrt(length(goodIndex));      
end
OKmeanList= meanPerColor;
OKerrList= errCO2;
OKmeanListPostCO2= meanPerColorPostCO2;
OKerrListPostCO2= errPostCO2;
OKmeanListAIR= meanPerColorAIR;
OKerrListAIR= errAIR;

pp2=subplot(2,2,2);
shadedErrorBar(x, OKmeanListAIR', OKerrListAIR','lineprops','-b')
%pp1.YLim= [0 max(max(meanList))+2000];
pp2.YLim= [0 150];
hold on;
shadedErrorBar(x, OKmeanList', OKerrList','lineprops','-r')
shadedErrorBar(x, OKmeanListPostCO2', OKerrListPostCO2','lineprops','-g')
hold off;
xlabel('Minute #');
%xlim([1, 120]);
ylabel(' % Counted ');
title(' % counts without 1 day' );
legend(lg, 'Location','eastoutside');

% 3.
for w=1:length(weirdIndex)
    sumPerColor=cumsum(p1(weirdIndex(w),1:duration)+p2(weirdIndex(w),1:duration), 2);
    sumPerColorPostCO2=cumsum(p1PostCO2(weirdIndex(w),1:duration)+p2PostCO2(weirdIndex(w),1:duration),2);
    sumPerColorAIR= cumsum(p1AIR(weirdIndex(w),1:duration)+p2AIR(weirdIndex(w),1:duration),2);
    pp3=subplot(2,2,w+2);
    plot(sumPerColorAIR);
    hold on
    plot(sumPerColor);
    plot(sumPerColorPostCO2);
    hold off
    xlabel('Minute #');
    ylabel(' Counts ');
    title(filesName(weirdIndex(w)));
    legend(lg, 'Location','eastoutside');
end    
% ==================================================
% ==================================================


% ==================================================   
% ==================================================
% Check the PostCO2 behavior of 5 hours exp vs 3 hours exp
% ==================================================
meanListPostCO2_5h= zeros(length(colorCues), duration);
errListPostCO2_5h= zeros(length(colorCues), duration);
meanListPostCO2_3h= zeros(length(colorCues), duration);
errListPostCO2_3h= zeros(length(colorCues), duration);
changeExp= 13; %index for 1st experiment of 3 hours duration
i=1;        
for c=colorCues
    indexes= find(strcmp(testedColorList, c));
    indexes5h= indexes(indexes < changeExp);
    indexes3h= indexes(indexes >= changeExp);
    %Sum the counts per position for each of the time groups
    sumPercolor5h=cumsum(p1PostCO2(indexes5h,1:duration)+p2PostCO2(indexes5h,1:duration), 2);
    sumPercolor3h=cumsum(p1PostCO2(indexes3h,1:duration)+p2PostCO2(indexes3h,1:duration),2);

    tempNorm= bsxfun(@rdivide,(sumPercolor5h*100),sumPercolor5h(:,end));
    tempNorm(isnan(tempNorm))=0;
    meanPerColorPostCO2_5h= mean(tempNorm);
    errPostCO2_5h= std(tempNorm)/sqrt(length(indexes5h));

    tempNorm= bsxfun(@rdivide,(sumPercolor3h*100),sumPercolor3h(:,end));
    tempNorm(isnan(tempNorm))=0;
    meanPerColorPostCO2_3h= mean(tempNorm);
    errPostCO2_3h= std(tempNorm)/sqrt(length(indexes3h));

    meanListPostCO2_5h(i,:)= meanPerColorPostCO2_5h;
    errListPostCO2_5h(i,:)= errPostCO2_5h;
    meanListPostCO2_3h(i,:)= meanPerColorPostCO2_3h;
    errListPostCO2_3h(i,:)= errPostCO2_3h;

    disp(indexes5h);
    disp(indexes3h);
    i=i+1;
end

%SHADED ERROR BAR
x= 1:duration; %[1:60; 1:60; 1:60; 1:60]
%Plot the counts withCO2 vs postCO2 each plot is for a value of gray
figure()
title('Counts in volume withCO2 vs postCO2 (per color)');
%set(gcf,'Title', title);
lg= [{'5h exp'}, {'3h exp'}];
for i= 1:4
    pp1=subplot(2,2,i);
    shadedErrorBar(x, meanListPostCO2_5h(i,:)', errListPostCO2_5h(i,:)','lineprops','-b')
    %pp1.YLim= [0 max(max(meanList))+2000];
    pp1.YLim= [0 150];
    hold on;
    shadedErrorBar(x, meanListPostCO2_3h(i,:)', errListPostCO2_3h(i,:)','lineprops','-g')
    hold off;
    xlabel('Minute #');
    %xlim([1, 120]);
    ylabel(' % Counted ');
    title(strcat(colorCues(i), {' '}, '(postCO2)'));
    legend(lg, 'Location','eastoutside');
end

% ==================================================
% (SAME WITH RAW DATA)
% Check the PostCO2 behavior of 5 hours exp vs 3 hours exp
% ==================================================   
changeExp= 13; %index for 1st experiment of 3 hours duration
i=1;
for c=colorCues
    indexes= find(strcmp(testedColorList, c));
    indexes5h= indexes(indexes < changeExp);
    indexes3h= indexes(indexes >= changeExp);
    lg5h= cell(1,length(indexes5h));
    lg3h= cell(1,length(indexes3h));
    %Sum the counts per position for each of the time groups
    sumPerColor5h=cumsum(p1PostCO2(indexes5h,1:duration)+p2PostCO2(indexes5h,1:duration), 2);
    sumPerColor3h=cumsum(p1PostCO2(indexes3h,1:duration)+p2PostCO2(indexes3h,1:duration),2);
    figure()

    for k= 1:2
        pp1=subplot(1,2,k);
        if k == 1
            plot(sumPerColor5h')
            %legend(indexes5h, 'Location','eastoutside');
            title(strcat(colorCues(i), {' '}, 'postCO2 (5h)'));
            lg5h= filesName(indexes5h);
            legend(lg5h, 'Location','eastoutside');

        else
            plot(sumPerColor3h')
            %legend(indexes3h, 'Location','eastoutside');
            title(strcat(colorCues(i), {' '}, 'postCO2 (3h)'));
            lg3h= filesName(indexes3h);
            legend(lg3h, 'Location','eastoutside');

        end 
        pp1.YLim= [0 max(max(sumPerColor5h(:,end)), max(sumPerColor3h(:,end)))];
        xlabel('Minute #');
        ylabel(' Counted (cumSum) ');
    end
    i= i+1;

end




%============================================================
% Working with time durations near a visual cue. 
% --> The insect behavior itis very similar to their behavior when 
%       counting frames where the
%       insect has been near the cue

% ========   
% * Do the dynamics change between experiments
%load the type of grays
colorCues= unique(testedColorList)';
% Initialize variable to obtain the mean and std error per colorCue  
i=1;
duration=60;
meanTimeListAIRp1= zeros(length(colorCues), duration);
errTimeListAIRp1= zeros(length(colorCues), duration);
meanTimeListCO2p1= zeros(length(colorCues), duration);
errTimeListCO2p1= zeros(length(colorCues), duration);
meanTimeListPostCO2p1= zeros(length(colorCues), duration);
errTimeListPostCO2p1= zeros(length(colorCues), duration);
meanTimeListAIRp2= zeros(length(colorCues), duration);
errTimeListAIRp2= zeros(length(colorCues), duration);
meanTimeListCO2p2= zeros(length(colorCues), duration);
errTimeListCO2p2= zeros(length(colorCues), duration);
meanTimeListPostCO2p2= zeros(length(colorCues), duration);
errTimeListPostCO2p2= zeros(length(colorCues), duration);

for c=colorCues
    indexes= find(strcmp(testedColorList, c));
    % Sum the times spend by some insects in the volume, is there any
    % difference with the counts plot?
    sumTimeInPosCO2p1= cumsum(t1CO2(indexes, 1:duration), 2);
    sumTimeInPosAIRp1= cumsum(t1AIR(indexes, 1:duration), 2);
    sumTimeInPosPostCO2p1= cumsum(t1PostCO2(indexes, 1:duration), 2);
    sumTimeInPosCO2p2= cumsum(t2CO2(indexes, 1:duration), 2);
    sumTimeInPosAIRp2= cumsum(t2AIR(indexes, 1:duration), 2);
    sumTimeInPosPostCO2p2= cumsum(t2PostCO2(indexes, 1:duration), 2);
    % obtain the total of time spent in both cues
    totalTimeCO2= sumTimeInPosCO2p1(:,end)+sumTimeInPosCO2p2(:,end);

    %Working with the time matrices for position 1
    tempNorm= bsxfun(@rdivide,(sumTimeInPosAIRp1*100),totalTimeCO2);
    tempNorm(isnan(tempNorm))=0;
    meanTimeAIRp1= mean(tempNorm);
    errTimeAIRp1= std(tempNorm)/sqrt(length(indexes));
    tempNorm= bsxfun(@rdivide,(sumTimeInPosCO2p1*100),totalTimeCO2);
    tempNorm(isnan(tempNorm))=0;
    meanTimeCO2p1= mean(tempNorm);
    errTimeCO2p1= std(tempNorm)/sqrt(length(indexes));
    tempNorm= bsxfun(@rdivide,(sumTimeInPosPostCO2p1*100),totalTimeCO2);
    tempNorm(isnan(tempNorm))=0;
    meanTimePostCO2p1= mean(tempNorm);
    errTimePostCO2p1= std(tempNorm)/sqrt(length(indexes));
    %Working with the time matrices for position 2
    tempNorm= bsxfun(@rdivide,(sumTimeInPosAIRp2*100),totalTimeCO2);
    tempNorm(isnan(tempNorm))=0;
    meanTimeAIRp2= mean(tempNorm);
    errTimeAIRp2= std(tempNorm)/sqrt(length(indexes));
    tempNorm= bsxfun(@rdivide,(sumTimeInPosCO2p2*100),totalTimeCO2);
    tempNorm(isnan(tempNorm))=0;
    meanTimeCO2p2= mean(tempNorm);
    errTimeCO2p2= std(tempNorm)/sqrt(length(indexes));
    tempNorm= bsxfun(@rdivide,(sumTimeInPosPostCO2p2*100),totalTimeCO2);
    tempNorm(isnan(tempNorm))=0;
    meanTimePostCO2p2= mean(tempNorm);
    errTimePostCO2p2= std(tempNorm)/sqrt(length(indexes));
    % Add the values calculated to its particular list
    meanTimeListCO2p1(i,:)= meanTimeCO2p1;
    errTimeListCO2p1(i,:)= errTimeCO2p1;
    meanTimeListAIRp1(i,:)= meanTimeAIRp1;
    errTimeListAIRp1(i,:)= errTimeAIRp1;
    meanTimeListPostCO2p1(i,:)= meanTimePostCO2p1;
    errTimeListPostCO2p1(i,:)= errTimePostCO2p1;
    meanTimeListCO2p2(i,:)= meanTimeCO2p2;
    errTimeListCO2p2(i,:)= errTimeCO2p2;
    meanTimeListAIRp2(i,:)= meanTimeAIRp2;
    errTimeListAIRp2(i,:)= errTimeAIRp2;
    meanTimeListPostCO2p2(i,:)= meanTimePostCO2p2;
    errTimeListPostCO2p2(i,:)= errTimePostCO2p2;
    i=i+1;
end

%USING SHADED ERROR BAR
x= 1:duration; %[1:60; 1:60; 1:60; 1:60]
%Plot the counts withCO2 vs postCO2 each plot is for a value of gray
figure(1);
title('Counts in volume withCO2 vs postCO2 (per color)');
%set(gcf,'Title', title);
lg= [{'Position 1'}, {'Position 2'}, {'Pos 1 + Pos 2'}];
for i= 1:4
    figure(1);
    pp1=subplot(2,2,i);
    % time near p1 and p2
    shadedErrorBar(x, meanTimeListAIR(i,:)', errTimeListAIR(i,:)','lineprops','-g')
    %pp1.YLim= [0 max(max(meanList))+2000];
    pp1.YLim= [0 100];
    hold on;
    shadedErrorBar(x, meanTimeListAIRp1(i,:)', errTimeListAIRp1(i,:)','lineprops','-b')
    shadedErrorBar(x, meanTimeListAIRp2(i,:)', errTimeListAIRp2(i,:)','lineprops','-r')

    hold off;
    xlabel('Minute #');
    %xlim([1, 120]);
    ylabel(' % time per position ');
    title(colorCues(i));
    legend(lg, 'Location','eastoutside');  
end

%total time spent in a cue per part of experiment
zAIR= sum(meanTimeListAIR,2);
zCO2= sum(meanTimeListCO2,2);
zPostCO2= sum(meanTimeListPostCO2, 2);
bar([zAIR, zCO2, zPostCO2])
xNames= [{'gray1'}, {'gray 2'}, {'gray 3'}, {'gray 4'}];
set(gca,'XtickLabel', xNames(1,:));



% ==================================================
% Estimate the time spent in each color cue (instead of position)
% ==================================================  
% Initialize matrices
i=1;
meanTimeListAIRcB= zeros(length(colorCues), duration);
errTimeListAIRcB= zeros(length(colorCues), duration);
meanTimeListAIRcT= zeros(length(colorCues), duration);
errTimeListAIRcT= zeros(length(colorCues), duration);

meanTimeListCO2cB= zeros(length(colorCues), duration);
errTimeListCO2cB= zeros(length(colorCues), duration);
meanTimeListCO2cT= zeros(length(colorCues), duration);
errTimeListCO2cT= zeros(length(colorCues), duration);

meanTimeListPostCO2cB= zeros(length(colorCues), duration);
errTimeListPostCO2cB= zeros(length(colorCues), duration);
meanTimeListPostCO2cT= zeros(length(colorCues), duration);
errTimeListPostCO2cT= zeros(length(colorCues), duration);

for c=colorCues
    % Find the experiment that used c as type of gray
    indexes= find(strcmp(testedColorList, c));
    % Find in which position the base and test color wehre placed
    p1BaseC= indexes(find(baseColorIndexList(indexes) == 1));
    p2BaseC= indexes(find(baseColorIndexList(indexes)== 2));

    %sum the times for each type ovf visual cue     
    sumTimeInBaseClrAIR= cumsum([t1AIR(p1BaseC, 1:duration); t2AIR(p2BaseC, 1:duration)], 2);
    sumTimeInTestClrAIR= cumsum([t1AIR(p2BaseC, 1:duration); t2AIR(p1BaseC, 1:duration)], 2);
    sumTimeInBaseClrCO2= cumsum([t1CO2(p1BaseC, 1:duration); t2CO2(p2BaseC, 1:duration)], 2);
    sumTimeInTestClrCO2= cumsum([t1CO2(p2BaseC, 1:duration); t2CO2(p1BaseC, 1:duration)], 2);
    sumTimeInBaseClrPostCO2= cumsum([t1PostCO2(p1BaseC, 1:duration); t2PostCO2(p2BaseC, 1:duration)], 2);
    sumTimeInTestClrPostCO2= cumsum([t1PostCO2(p2BaseC, 1:duration); t2PostCO2(p1BaseC, 1:duration)], 2);

    % obtain the total of time spent in both cues
    totalTimeCO2= sumTimeInBaseClrCO2(:,end) +sumTimeInTestClrCO2(:,end);

    %Normalize the amount of time spend near each color cue
    tempNorm= bsxfun(@rdivide,(sumTimeInBaseClrAIR*100),totalTimeCO2);
    tempNorm(isnan(tempNorm))=0;
    meanTimeBaseClrAIR= mean(tempNorm);
    errTimeBaseClrAIR= std(tempNorm)/sqrt(length(indexes));
    tempNorm= bsxfun(@rdivide,(sumTimeInTestClrAIR*100),totalTimeCO2);
    tempNorm(isnan(tempNorm))=0;
    meanTimeTestClrAIR= mean(tempNorm);
    errTimeTestClrAIR= std(tempNorm)/sqrt(length(indexes));     

    tempNorm= bsxfun(@rdivide,(sumTimeInBaseClrCO2*100),totalTimeCO2);
    tempNorm(isnan(tempNorm))=0;
    meanTimeBaseClrCO2= mean(tempNorm);
    errTimeBaseClrCO2= std(tempNorm)/sqrt(length(indexes));
    tempNorm= bsxfun(@rdivide,(sumTimeInTestClrCO2*100),totalTimeCO2);
    tempNorm(isnan(tempNorm))=0;
    meanTimeTestClrCO2= mean(tempNorm);
    errTimeTestClrCO2= std(tempNorm)/sqrt(length(indexes));

    tempNorm= bsxfun(@rdivide,(sumTimeInBaseClrPostCO2*100),totalTimeCO2);
    tempNorm(isnan(tempNorm))=0;
    meanTimeBaseClrPostCO2= mean(tempNorm);
    errTimeBaseClrPostCO2= std(tempNorm)/sqrt(length(indexes));
    tempNorm= bsxfun(@rdivide,(sumTimeInTestClrPostCO2*100),totalTimeCO2);
    tempNorm(isnan(tempNorm))=0;
    meanTimeTestClrPostCO2= mean(tempNorm);
    errTimeTestClrPostCO2= std(tempNorm)/sqrt(length(indexes));

    % Add values to their corresponding list to plot
    meanTimeListAIRcB(i,:)= meanTimeBaseClrAIR;
    errTimeListAIRcB(i,:)= errTimeBaseClrAIR;
    meanTimeListAIRcT(i,:)= meanTimeTestClrAIR;
    errTimeListAIRcT(i,:)= errTimeTestClrAIR;

    meanTimeListCO2cB(i,:)= meanTimeBaseClrCO2;
    errTimeListCO2cB(i,:)= errTimeBaseClrCO2;
    meanTimeListCO2cT(i,:)= meanTimeTestClrCO2;
    errTimeListCO2cT(i,:)= errTimeTestClrCO2;

    meanTimeListPostCO2cB(i,:)= meanTimeBaseClrPostCO2;
    errTimeListPostCO2cB(i,:)= errTimeBaseClrPostCO2;
    meanTimeListPostCO2cT(i,:)= meanTimeTestClrPostCO2;
    errTimeListPostCO2cT(i,:)= errTimeTestClrPostCO2;

    i= i+1;
end
clear meanTimeBaseClrAIR errTimeBaseClrAIR meanTimeTestClrAIR errTimeTestClrAIR
clear meanTimeBaseClrCO2 errTimeBaseClrCO2 meanTimeTestClrtCO2 errTimeTestClrCO2
clear meanTimeBaseClrPostCO2 errTimeBaseClrPostCO2 meanTimeTestClrPostCO2 errTimeTestClrPostCO2



%USING SHADED ERROR BAR
x= 1:duration; %[1:60; 1:60; 1:60; 1:60]
%Plot the counts withCO2 vs postCO2 each plot is for a value of gray
figure(1);
%set(gcf,'Title', title);
lg= [{'Base color'}, {'Test color'}];
for i= 1:4
    figure(1);
    pp1=subplot(2,2,i);
    % time near p1 and p2
    shadedErrorBar(x, meanTimeListCO2cB(i,:)', errTimeListCO2cB(i,:)','lineprops','-g')
    %pp1.YLim= [0 max(max(meanList))+2000];
    pp1.YLim= [0 100];
    hold on;
    shadedErrorBar(x, meanTimeListCO2cT(i,:)', errTimeListCO2cT(i,:)','lineprops','-b')

    hold off;
    xlabel('Minute #');
    %xlim([1, 120]);
    ylabel(' % time per position ');
    title(colorCues(i));
    legend(lg, 'Location','eastoutside');  
end
%Data from AIR
figure(2)
for i= 1:4
    figure(2);
    pp1=subplot(2,2,i);
    % time near p1 and p2
    shadedErrorBar(x, meanTimeListAIRcB(i,:)', errTimeListAIRcB(i,:)','lineprops','-g')
    %pp1.YLim= [0 max(max(meanList))+2000];
    pp1.YLim= [0 100];
    hold on;
    shadedErrorBar(x, meanTimeListAIRcT(i,:)', errTimeListAIRcT(i,:)','lineprops','-b')
    hold off;
    xlabel('Minute #');
    %xlim([1, 120]);
    ylabel(' % time per position ');
    title(colorCues(i));
    legend(lg, 'Location','eastoutside');  
end
%Data from POSTCO2
figure(3)
for i= 1:4
    figure(3);
    pp1=subplot(2,2,i);
    % time near p1 and p2
    shadedErrorBar(x, meanTimeListPostCO2cB(i,:)', errTimeListPostCO2cB(i,:)','lineprops','-g')
    %pp1.YLim= [0 max(max(meanList))+2000];
    pp1.YLim= [0 100];
    hold on;
    shadedErrorBar(x, meanTimeListPostCO2cT(i,:)', errTimeListPostCO2cT(i,:)','lineprops','-b')
    hold off;
    xlabel('Minute #');
    %xlim([1, 120]);
    ylabel(' % time per position ');
    title(colorCues(i));
    legend(lg, 'Location','eastoutside');  
end
% ==================================================
% ==================================================



% ==================================================   
% PLOT SINGLE TRAJ OF THE INSECT THAT STOOD THE MOST IN A CUE
% You can load the information from a particular file and plot the
% trajectory with the largest time in a cue
% ==================================================
filesPath= strcat(outputPath, outputFolder);
numGrps=120; %=120 =60 sec per group || =4 =1800 sec per grp
% Choose que file related to the dataset do you want to open
fileIndex= 5;
fileName= strcat(dataset(fileIndex).fileName(1:15),'_countsInsideCueVol_CO2.xlsx');

%Matrices with the different ts values [(t1-t0), ..., (tn- t(n-1))] 
tDiff1 = zeros(length(filesList),numGrps); 
tDiff2 = zeros(length(filesList),numGrps);

% load the counts from the .xlsx file
fileName= filesList(fileIndex).name;
disp(strcat(' - Working with file: ', {' '}, fileName));
dataFromExcel = table2array(readtable(strcat(filesPath, fileName)));
initialTS= dataFromExcel(1,4);    
finalTS= dataFromExcel(1,5);

%sort matrix in fct of the cue-id-timestamp
sortedData= sortrows(dataFromExcel(:,1:3));
%find the first appearence of the position 2 
split= find(sortedData(:,1) ==2,1);


%counts the number of repetiions for each ID in each cue
%Cue in Pos 1
[IDsCtsP1,IDsGrpP1] = groupcounts(sortedData(1:(split-1),2));
%Cue in Pos 2
[IDsCtsP2,IDsGrpP2] = groupcounts(sortedData(split:end,2));
%Find the estimation of the total duration spent by all insects in each cue
timeInPos1= sum(IDsCtsP1)/fps
timeInPos2= sum(IDsCtsP2)/fps 
% preference index towards POS 1
prefIndexPos=(timeInPos1 - timeInPos2)/ (timeInPos1+timeInPos2)


%Find Trajectories IDs with biger trajectories (at least half the
%biggest trajectory)

% Find the longest trajectory
ix= find(IDsCtsP1 == max(IDsCtsP1));
idToLook= IDsGrpP1(ix)';
ix= find(IDsCtsP2 == max(IDsCtsP2));
idToLook= horzcat(idToLook, IDsGrpP2(ix)');

%Select wichi ID do you want to plot
indexesID= find(dataset(fileIndex).attr_id(:) == idToLook(1));
%load experiment cues positions
cuesSetup= dataset(fileIndex).expCues;
% parameters for size of the volume around the cue and wind tunnel dimensions
radius= 0.07;
h= 0.04;
lim_x= 0.9144;
subplot(2,1,1);
plt1=plot(dataset(fileIndex).attr_x(indexesID), dataset(fileIndex).attr_y(indexesID), 'b');
hold on
% add starting and ending points
plot(dataset(fileIndex).attr_x(indexesID(1)), dataset(fileIndex).attr_y(indexesID(1)),'b', 'Marker', '>');
plot(dataset(fileIndex).attr_x(indexesID(end)), dataset(fileIndex).attr_y(indexesID(end)), 'xb');
%build the volume
x= -lim_x  + cell2mat(cuesSetup(2,2));
y= cell2mat(cuesSetup(2,3));
x1= x - radius;
x2= x + radius;
y1= y - radius;
y2= y + radius;
xValues= [x1, x2, x2, x1, x1];
yValues= [y1, y1, y2, y2, y1];
plot(xValues,yValues, ':r');
x= -lim_x  + cell2mat(cuesSetup(3,2));
y= cell2mat(cuesSetup(3,3));
x1= x - radius;
x2= x + radius;
y1= y - radius;
y2= y + radius;
xValues= [x1, x2, x2, x1, x1];
yValues= [y1, y1, y2, y2, y1];
plot(xValues,yValues, ':k');
plt1= add_visual_cues_to_plot(plt1, cuesSetup(:,[1 2 3]));
hold off
xlim([-0.91 0.91]);
ylim([-0.30 0.30]);
ylabel(' Y ');
xlabel(' X '); 
title(strcat('CO2 Trajectory: ',{' '}, num2str(idToLook(1)), {' '}, ' -Exp:', {' '}, dataset(fileIndex).fileName(1:8), '-', dataset(fileIndex).fileName(10:15)));   
subplot(2,1,2);
plt2=plot(dataset(fileIndex).attr_x(indexesID), dataset(fileIndex).attr_z(indexesID), 'b');
hold on
% add starting and ending points
plot(dataset(fileIndex).attr_x(indexesID(1)), dataset(fileIndex).attr_z(indexesID(1)),'b', 'Marker', '>');
plot(dataset(fileIndex).attr_x(indexesID(end)), dataset(fileIndex).attr_z(indexesID(end)), 'xb');
%build the volume
x= -lim_x  + cell2mat(cuesSetup(3,2));
z= 0;
x1= x - radius;
x2= x + radius;
z1= z ;
z2= z + h;
xValues= [x1, x2, x2, x1, x1];
zValues= [z1, z1, z2, z2, z1];
plot(xValues,zValues, ':k');
plt2= add_visual_cues_to_plot(plt2, cuesSetup(:,[1 2 4]));
hold off
xlim([-0.91 0.91]);
ylim([0.00 0.6]);    
ylabel(' Z ');
xlabel(' X ');
% =============================================  
  
    
