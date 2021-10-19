
%PI grouped by mutant type
% With std error bar
for expIndex=1:length(dataset)
    typeGenderList{expIndex}= strcat(dataset(expIndex).type,'-',dataset(expIndex).gender);
end
uniqueTG= unique(typeGenderList);
tgIndex=1;
for tgValue=uniqueTG
    %disp(tgValue);
    i= find(strcmp(typeGenderList, tgValue));
    meanPIperTG(tgIndex)= mean(piList(i));
    %stdError/number of experiments for this type/gender
    errList(tgIndex)= std(piList(i))/length(i);
    tgIndex= tgIndex+1;
end

% Sort the Data & Rearrange Labels
[sortedPIs, newIndexes] = sort(meanPIperTG); % sorts in *ascending* order

sortedLabels = uniqueTG(newIndexes); 
sortedErr= errList(newIndexes);
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


%===========================================================
 % as scatttered points
for expIndex=1:length(dataset)
    typeGenderList{expIndex}= strcat(dataset(expIndex).type,'-',dataset(expIndex).gender);
end
uniqueTG= unique(typeGenderList);
tgIndex=1;
prefIndexesTG=[];
xAxis= [];
for tgValue=uniqueTG
    %disp(tgValue);
    i= find(strcmp(typeGenderList, tgValue));
    prefIndexesTG=vertcat(prefIndexesTG, piList(i));
    meanPIperTG(tgIndex)= mean(piList(i));

    xAxis= vertcat(xAxis, zeros(length(i),1)+tgIndex);
    tgIndex= tgIndex+1;
end

axes('Xlim', [0, length(piList)+1], 'XTick', 0:20:length(piList)+1);
scatter(xAxis, prefIndexesTG, 'jitter', 'on');
hold on
plot(meanPIperTG, '+r');
hold off




% tgIndex=1;
% prefIndexesTG=[];
% xAxis= [];
% 
% tgValue='l4-f';
% i= find(strcmp(typeGenderList, tgValue));
% prefIndexesTG=vertcat(prefIndexesTG, piList(i));
% xAxis= vertcat(xAxis, zeros(length(i),1)+tgIndex);
% tgIndex= tgIndex+1;



% Sort the Data & Rearrange Labels
[sortedPIs, newIndexes] = sort(meanPIperTG); % sorts in *ascending* order

sortedLabels = uniqueTG(newIndexes); 
sortedErr= errList(newIndexes);
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
 
