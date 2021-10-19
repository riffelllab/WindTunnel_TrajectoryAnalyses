% Plot the number of insect detected near each of the visual cues
% Arguments:
%   - data: matrix withe the number of insects counted in each volume per
%           visual cue (columns) and per experiments (row)                      
%	        [InsectCountedInCue1, ..., InsectCountedInCueN]
%   - subFolder: Indicator of the experiment being currently analyzed
%   - checkedType: Type of mosquitoes being currently analyzed (wt,l...)
%   - checkedGender: Gender of the mosquitoes currently being analyzed
function plot_insect_near_position(data, checkedType, checkedGender)
    % Calculate the total amount of trajectories in the FULL dataset
    totalTraj= sum(sum(data));
    
    % Plot the Mean counts for the each of XY positions of the V.Cues for all '4_colors' experiments
    figure()

    if length(data(1,:)) == 4
        disp('4 positions');
        bar([sum(data(:,1))/totalTraj; sum(data(:,2))/totalTraj; sum(data(:,3))/totalTraj; sum(data(:,4))/totalTraj]);
        xV= cell(1,4);
        xV{1}='Pos 1 (0.45, 0.10)'; xV{2}='Pos 2 (0.65, 0.10)'; xV{3}='Pos 3 (0.45, -0.10)';xV{4}='Pos 1 (0.65, -0.10)';
    else 
        bar([sum(data(:,1))/totalTraj; sum(data(:,2))/totalTraj]);
        xV= cell(1,2);
        xV{1}='Pos 1'; xV{2}='Pos 2';
        %xV{1}='Pos 2 (0.67, 0.095)'; xV{2}='Pos 4 (0.67, -0.095)';
    end
    set(gca, 'XtickLabel', xV);
    ylabel('% per position');

    if contains(version, 'R2019')
        xtickangle(45)    
    end
    title(strcat('Counts for: ',{' '}, checkedType, {' - '}, checkedGender));

    %print(gcf,strcat('count_near_position.tif'),'-dtif','-r300');        % *// 300 dpi

end