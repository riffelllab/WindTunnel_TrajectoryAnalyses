
% Plots counts per position and per experiment. Data contains the counts 
% for all the experiments (one experiment day per row).
% WARNING: Colors assigned to position must be added a posteriori (via Paint)
% Arguments:
%   - data: [expIndex countsForPos_1 countsForPos_2 countsForPos_3 countsForPos_4]
%   - xNames: List with the files names used for data
%   - subFolder: Indicator of the experiment being currently analyzed
function plot_insect_near_position_per_exp_over_time(data, xNames, subFolder)

figure()
bp= bar(data(:,:))


set(gca,'XtickLabel', xNames);
if contains(version, 'R2019')
    xtickangle(45)    
end
%create the legend in function of the type of experiment
lg= cell(1,4);
if contains(subFolder, '4_grays')
    lg{1}='POS-1'; lg{2}='POS-2'; lg{3}='POS-3';lg{4}='POS-4';    
elseif contains(subFolder, '4_reds')
    lg{1}='R-T1 '; lg{2}='RW-HUE'; lg{3}='R-HUE';lg{4}='RO-S1';
elseif contains(subFolder, '4_colors')
    lg{1}='WHITE'; lg{2}='GREEN'; lg{3}='RED';lg{4}='BLUE';
elseif (contains(subFolder, '2_blacks') || contains(subFolder, 'white_red') || contains(subFolder, 'white_gray'))
    lg= cell(1,2);
    lg{1}='POS-1'; lg{2}='POS-2';
else
    lg{1}='POS-1'; lg{2}='POS-2'; lg{3}='POS-3';lg{4}='POS-4';
end
legend(bp,lg, 'Location', 'northwest');
title(' Counts in volume near position per experiment');

end