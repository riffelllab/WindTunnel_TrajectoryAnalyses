
%Function to add the odor and visual cues used in the experiment to the plot.
%Arguments:
%   - figPlot:      Current plot where we want to add the odor and visual cues
%   - cluesSetup:   Cell matrix with the odor/visual cues used and their XYZ
%                   position
%Returns:
%   - figPlot:      Current plot with the cues added
function [figPlot]= add_visual_cues_to_plot(figPlot, cluesSetup)

    hold on;
    %Load the clues used (if any)
    if ~isempty(cluesSetup)
        %Load the XY position for ODOR
        x= -0.9 + cell2mat(cluesSetup(1,2));
        y= cell2mat(cluesSetup(1,3));

        % Plot a start where it is the Odor source
        plot(x,y, 'p', 'markerfacecolor','magenta', 'MarkerSize', 8);

        for i= 2:length(cluesSetup(:,1))
            % Load the XY position for each visual clue
            x=-0.9 + cell2mat(cluesSetup(i,2));
            y= cell2mat(cluesSetup(i,3));
            % Select the correct Matlab color code for the coler used as v. cue
            vCueColor= convert_color(char(cluesSetup(i,1)));
            % Plot the visual cue
            plot(x, y, 'o', 'MarkerSize', 10, 'markerfacecolor', vCueColor);
        end     
    end
    hold off;
end