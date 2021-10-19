% Plot the number of insect detected near each of the visual cues

function plot_insect_counted_in_vol(data, subFolder, checkedType, checkedGender)
    % Calculate the total amount of trajectories in the FULL dataset
    totalTraj= sum(sum(data(:,2:end)));
    
    % Plot the Mean count for the XY positions of the V.Cues for all '4_colors' experiments
    figure()
    
    if  ~isempty(strfind(subFolder, 'black_white'))
        % If subFolder contains 'black_white' use these values for the X axis
        xV= cell(1,6);
        xV{1}='Black Clue (m0 - f)'; xV{2}='White Clue (m0 - f)'; xV{3}='Black Clue (m1 - f)'; xV{4}='White Clue (m1 - f)';xV{5}='Black Clue (wt - f)'; xV{6}='White Clue (wt - f)';

    elseif ~isempty(strfind(subFolder, '4_colors'))
        % If subFolder contains '4_colors' use these values for the X axis
        disp('Plotting counts for 4 colors')
        % Add all the counts for each color trhough all the experiments
        % Add all the counts for the white cue
        whiteCounts= data(1,2)+ data(2, 4)+ data(3,3)+data(4,5);
        % Add all the counts for the green cue
        greenCounts= data(1,3)+ data(2, 5)+ data(3,2)+data(4,4);
        % Add all the counts for the red cue
        redCounts= data(1,4)+ data(2, 2)+ data(3,5)+data(4,3);
        % Add all the counts for the blue cue
        blueCounts= data(1,5)+ data(2, 3)+ data(3,4)+data(4,2);

        % Plot the Mean count for the V.Cues (regardless their position) for all '4_colors' experiments
        %bar([mean(countListCO2(1:2, 2)); mean(countListCO2(1:2, 3)); mean(countListCO2(4:5, 2)); mean(countListCO2(4:5, 3)); mean(countListCO2(6, 2)); mean(countListCO2(6, 3))])
        bar([whiteCounts/totalTraj; greenCounts/totalTraj; redCounts/totalTraj; blueCounts/totalTraj]);
        xV= cell(1,4);
        xV{1}='White cue'; xV{2}='Green cue'; xV{3}='Red cue';xV{4}='Blue cue';

    elseif ~isempty(strfind(subFolder, '4_reds'))
        % If subFolder contains '4_colors' use these values for the X axis
        disp('plotting counts for 4 types of red')
        % Add all the counts for each red (R-T1, RW-HUE, R-HUE, RO-S1) trhough all the experiment
        RT1counts= sum(data(1, 2)+ data(2, 3)+data(3, 5)+ data(4, 4));
        RWHUEcounts= sum(data(1, 3)+ data(2, 5)+data(3, 4)+ data(4, 2));
        RHUEcounts= sum(data(1, 4)+ data(2, 2)+data(3, 3)+ data(4, 5));
        ROS1counts= sum(data(1, 5)+ data(2, 4)+data(3, 2)+ data(4, 3));

        bar([RT1counts/totalTraj; RWHUEcounts/totalTraj; RHUEcounts/totalTraj; ROS1counts/totalTraj]);
        xV= cell(1,4);
        xV{1}='R-T1 '; xV{2}='RW-HUE'; xV{3}='R-HUE';xV{4}='RO-S1';
    elseif ~isempty(strfind(subFolder, 'HistamineFed'))
        % If subFolder contains '4_colors' use these values for the X axis
        bar([data(1,2)/totalTraj; data(1,3)/totalTraj; data(1,5)/totalTraj; data(1,5)/totalTraj]);
        xV= cell(1,4);
        xV{1}=' GRAY-1 '; xV{2}='WHITE'; xV{3}='BLACK';xV{4}='GRAY-2';
    end;
    set(gca, 'XtickLabel', xV);
    title(strcat('Counts for: ',{' '}, checkedType, {' - '}, checkedGender));

end