%Function to plot a given insect trajectory in a 2D plot. It generates 2
%subplots: One showing the current  trajectory in the (X,Y) plane and
%another showing the current trajectory in the (X,Z) plane
%Arguments: 
%   - id: ID value for the insect trajectory to plot
%   - objXYZ: trajectory of the insect (XYZ position over time)
%   - color: Color tro use in the insect trajectory plot

function [fig]= plot_trajectory_2D_v3(objID, objXYZ, color, cluesSetup, duration, init)

    fig= figure(gcf);    
    subplot(2,1,1);
    %Plot the (X,Y) view of the insect trajectory
    plot(objXYZ(:,1), objXYZ(:,2), 'Color',color)
    hold on;
    %Mark the initial position for current objID
    plot(objXYZ(1,1),objXYZ(1,2),'-<', 'Color', color);
    l1= ''%strcat('Starting time: ',num2str(init), 's.');
    l2= ''%strcat(' - Duration: ',num2str(duration), 's.');
    text= sprintf(strcat('(X,Y) view of the test section. Insect ID: ',num2str(objID),'\n',l1,l2));
    %str = {l1,l2};
    %dim = [0, 0.5, 0, 0];
    %annotation('textbox', dim,'String',str,'FitBoxToText','on');
    %Mark the last position for current objID
    plot(objXYZ(end,1),objXYZ(end,2),'-x', 'Color',color);
    % Load the clues used (if any)
    if isempty(cluesSetup) == 0
        %Load the XY position for ODOR
        x= -0.9 + cell2mat(cluesSetup(1,2));
        y= cell2mat(cluesSetup(1,3));
        plot(x,y, 'p','markersize', 10, 'markerfacecolor','magenta');
        for i= 2:length(cluesSetup(:,1))
            %Load the XY position for each visual clue
            x= -0.9 + cell2mat(cluesSetup(i,2));
            y= cell2mat(cluesSetup(i,3));
            plot(x, y, 'o', 'markersize', 10, 'markerfacecolor', char(cluesSetup(i,1)));
        end;
    end;
    hold off;
    %title(strcat('(X,Y) view of the test section. Insect ID: ',num2str(objID)));
    title(text);
    xlabel('X axis');
    ylabel('Y axis');
    xlim([-0.91 0.91]);
    ylim([-0.30 0.30]);
    %set(gca, 'YDir','Reverse');
    
    subplot(2,1,2);
    %Plot the (X,Z) view of the insect trajectory
    plot(objXYZ(:,1), objXYZ(:,3), 'Color',color)
    hold on;
    %Mark the initial position for current objID
    plot(objXYZ(1,1),objXYZ(1,3),'-<', 'Color', color);    
    %Mark the last position for current objID
    plot(objXYZ(end,1),objXYZ(end,3),'-x', 'Color',color);
    %Load the clues used (if any)
    if isempty(cluesSetup) == 0
        %Load the XY position for ODOR
        x= -0.9 + cell2mat(cluesSetup(1,2));
        z= cell2mat(cluesSetup(1,4));
        plot(x,z, 'p','markersize', 10, 'markerfacecolor','magenta');
        for i= 2:length(cluesSetup(:,1))
            %Load the XY position for each visual clue
            x= -0.9 + cell2mat(cluesSetup(i,2));
            z= cell2mat(cluesSetup(i,4));
            plot(x, z, 'o', 'markersize', 10, 'markerfacecolor', char(cluesSetup(i,1)));
        end;
    end;
    title(strcat('(X,Z) view of the test section. Insect ID: ',num2str(objID)));
    %title(text);
    xlabel('X axis');
    ylabel('Z axis');
    xlim([-0.91 0.91]);
    ylim([0.0 0.60]);
    %ylim([0 0.60]);
    hold off;
    
    %Save the plot as a TIFF image
    print(gcf,strcat('traj_insect_',num2str(objID), '.tiff'),'-dtiff','-r300');        % *// 300 dpi

    %disp(strcat(' * objID over 3 seconds: ',num2str(objID)));
    
    
end