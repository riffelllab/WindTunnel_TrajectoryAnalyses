% Plot the number of insect detected near each of the visual cues
% Arguments:
%   - data: matrix withe the number of insects counted in each volume per
%           visual cue (columns) and per experiments (row)                      
%           [InsectCountedInCue1, ..., InsectCountedInCueN]
%   - subFolder: Indicator of the experiment being currently analyzed
%   - checkedType: Type of mosquitoes being currently analyzed (wt,l...)
%   - checkedGender: Gender of the mosquitoes currently being analyzed

function plot_insect_near_vCue(data, subFolder, checkedType, checkedGender)
    % Calculate the total amount of trajectories in the FULL dataset
    totalTraj= sum(sum(data(:,:)));
    
    % Plot the Mean count for the XY positions of the V.Cues for all '4_colors' experiments
    figure()
    
    if  contains(subFolder, 'mutants')
        WHITEcounts= data(1,1)+ data(2,2)+ data(3,2)+ data(4,1)+ data(5,1)+ data(6,2)+ data(7,1)+ data(8,1)+ data(9,2)+ ...
            data(10,1)+ data(11,2)+ data(12,2)+ data(13,1)+ data(14,2)+ data(15,2)+ data(16,1)+ data(17,1)+ data(18,1)+ data(19,2);        
        BLACKCounts= data(1,2)+ data(2,1)+ data(3,1)+ data(4,2)+ data(5,2)+ data(6,1)+ data(7,2)+ data(8,2)+ data(9,1)+ ...
            data(10,2)+ data(11,1)+ data(12,1)+ data(13,2)+ data(14,1)+ data(15,1)+ data(16,2)+ data(17,2)+ data(18,2)+ data(19,1);    
        
        bar([WHITEcounts/totalTraj; BLACKCounts/totalTraj]);
        xV= cell(1,2);
        xV{1}=' WHITE (mut)'; xV{2}='BLACK (mut)';    
        
    elseif contains(subFolder, 'black_white')
        % If subFolder contains 'black_white' use these values for the X axis
        xV= cell(1,6);
        xV{1}='Black Clue (m0 - f)'; xV{2}='White Clue (m0 - f)'; xV{3}='Black Clue (m1 - f)'; xV{4}='White Clue (m1 - f)';xV{5}='Black Clue (wt - f)'; xV{6}='White Clue (wt - f)';

    elseif contains(subFolder, '4_colors')
        % If subFolder contains '4_colors' use these values for the X axis
        disp('Plotting counts for 4 colors')
        % Add all the counts for each color trhough all the experiments
        % Add all the counts for the white cue
        whiteCounts= data(1,4)+ data(2, 3)+ data(3,1)+data(4,4)+data(5,2);
        % Add all the counts for the green cue
        greenCounts= data(1,3)+ data(2, 1)+ data(3,4)+data(4,2)+data(5,3);
        % Add all the counts for the red cue
        redCounts= data(1,2)+ data(2, 4)+ data(3,3)+data(4,1)+data(5,1);
        % Add all the counts for the blue cue
        blueCounts= data(1,1)+ data(2,2)+ data(3,2)+data(4,3)+data(5,4);

        % Plot the Mean count for the V.Cues (regardless their position) for all '4_colors' experiments
        %bar([mean(countListCO2(1:2, 2)); mean(countListCO2(1:2, 3)); mean(countListCO2(4:5, 2)); mean(countListCO2(4:5, 3)); mean(countListCO2(6, 2)); mean(countListCO2(6, 3))])
        bar([whiteCounts/totalTraj; greenCounts/totalTraj; redCounts/totalTraj; blueCounts/totalTraj]);
        xV= cell(1,4);
        xV{1}='WHITE'; xV{2}='GREEN'; xV{3}='RED';xV{4}='BLUE';

    elseif contains(subFolder, '4_reds')
        % If subFolder contains '4_colors' use these values for the X axis
        disp('plotting counts for 4 types of red')
        % Add all the counts for each red (R-T1, RW-HUE, R-HUE, RO-S1) trhough all the experiment
        RT1counts= sum(data(1,1)+ data(2,2)+data(3,4)+ data(4,3));
        RWHUEcounts= sum(data(1,2)+ data(2,4)+data(3,3)+ data(4,1));
        RHUEcounts= sum(data(1,3)+ data(2,1)+data(3,2)+ data(4,4));
        ROS1counts= sum(data(1,4)+ data(2,3)+data(3,1)+ data(4,2));

        bar([RT1counts/totalTraj; RWHUEcounts/totalTraj; RHUEcounts/totalTraj; ROS1counts/totalTraj]);
        xV= cell(1,4);
        xV{1}='R-T1 '; xV{2}='RW-HUE'; xV{3}='R-HUE';xV{4}='RO-S1';
        
    elseif contains(subFolder, '4_grays')
        % If subFolder contains '4_grays' use these values for the X axis
        GRAY1counts= sum(data(1,2)+ data(2,3)+data(3,1)+ data(4,4));
        GRAY2counts= sum(data(1,1)+ data(2,4)+data(3,2)+ data(4,3));
        GRAY3counts= sum(data(1,4)+ data(2,2)+data(3,3)+ data(4,1));
        GRAY4counts= sum(data(1,3)+ data(2,1)+data(3,4)+ data(4,2));
        
        bar([GRAY1counts/totalTraj; GRAY2counts/totalTraj; GRAY3counts/totalTraj; GRAY4counts/totalTraj]);
        xV= cell(1,4);
        xV{1}=' GRAY-1 '; xV{2}='GRAY-2'; xV{3}='GRAY-3';xV{4}='GRAY-4';
    
    elseif contains(subFolder, 'white_red')
        WHITEcounts= data(1,1)+ data(2,2)+ data(3,1)+ data(4,2)+ data(5,1);
        REDcounts= data(1,2)+ data(2,1) + data(3,2)+ data(4,1)+ data(5,2);
        
        bar([WHITEcounts/totalTraj; REDcounts/totalTraj]);
        xV= cell(1,2);
        xV{1}=' WHITE '; xV{2}='RED';     
    
    elseif contains(subFolder, 'white_gray')
        WHITEcounts= (data(1,2)+data(2,1)+data(3,1)+data(4,2)+data(5,2)+data(6,2)+data(8,1)+data(9,2)+data(10,1)+data(11,1)+data(12,1) ...
                      +data(13,1)+data(14,2)+data(15,2)+data(16,1)+data(17,1)+data(18,2)+data(19,1)+data(20,2));
        GRAY1counts= data(1,1)+data(5,1)+data(12,2)+data(13,2)+ data(14,1);
        GRAY2counts= data(6,1)+data(7,1)+data(9,1)+data(19,2)+data(20,1);
        GRAY3counts= data(2,2)+data(4,1)+data(10,2)+data(15,1)+data(16,2);
        GRAY4counts= data(3,2)+data(8,2)+data(11,2)+data(17,2)+data(18,1);
        totalGrayX= GRAY1counts+GRAY2counts+GRAY3counts+GRAY4counts;
        bar([GRAY1counts/totalTraj; GRAY2counts/totalTraj; GRAY3counts/totalTraj; GRAY4counts/totalTraj; WHITEcounts/totalTraj; totalGrayX/totalTraj ]);
        xV= cell(1,5);
        xV{1}='GRAY1'; xV{2}='GRAY2'; xV{3}='GRAY3'; xV{4}='GRAY4'; xV{5}='WHITE'; xV{6}='ALL GRAYS';
      
    elseif contains(subFolder, 'black_red')
        BLACKcounts= data(1,1)+ data(2,1)+ data(3,2)+ data(4,1)+ data(5,2)+ data(6,2)+ data(7,2);
        REDcounts= data(1,2)+ data(2,2) + data(3,1)+ data(4,2)+ data(5,1)+ data(6,1)+ data(7,1);
        
        bar([BLACKcounts/totalTraj; REDcounts/totalTraj]);
        xV= cell(1,2);
        xV{1}=' BLACK '; xV{2}='RED';   
    elseif contains(subFolder, 'gray_red')
        GRAYcounts= data(1,2)+ data(2,1)+ data(3,1)+ data(4,2);
        REDcounts= data(1,1)+ data(2,2) + data(3,2)+ data(4,1);
        
        bar([GRAYcounts/totalTraj; REDcounts/totalTraj]);
        xV= cell(1,2);
        xV{1}=' GRAY '; xV{2}='RED';       
    elseif contains(subFolder, 'Histamine')
        disp('    * PLOTTING DATA FOR 1 EXPERIMENT ONLY!!!');
        bar(data)
        xV= cell(1,2);
        xV{1}=' BLACK '; xV{2}=' WHITE ';  
    else
        disp('   * ERROR! Experiment plot not implemented yet')
    end
    set(gca, 'XtickLabel', xV);
    ylabel('% per visual cue');
    title(strcat('Counts for: ',{' '}, checkedType, {' - '}, checkedGender));

end