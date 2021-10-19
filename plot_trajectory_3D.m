
%Function to plot a given insect trajectory
%Arguments: 
%   - id: ID value for the insect trajectory to plot
%   - objXYZ: trajectory of the insect (XYZ position over time)
%   - testSectionVol: Plot figure containing the test section volumen. The trajectory of the insect will be plotted in this figure 
%   - color: Color tro use in the insect trajectory plot

function plot_trajectory_3D(id, objXYZ, testSectionVol, color, timeTreshold)

    %Create the test section in the plot to add the mosquito paths
    %testSectionVol= load_test_section_volumen();
    
    figure(testSectionVol);
    hold on
    %Plot the XYZ values of the current objID
    plot3(objXYZ(:,1),objXYZ(:,2),objXYZ(:,3), 'Color',color)
            
    %Mark the initial position for current objID
    plot3(objXYZ(1,1),objXYZ(1,2),objXYZ(1,3),'-o', 'Color', color);
        
    %Mark the last position for current objID
    plot3(objXYZ(end,1),objXYZ(end,2),objXYZ(end,3),'-x', 'Color',color);
    
    hold off
  %  grid on;
    
    title(strcat('Path tracked with flydra for insects with duration >  ',num2str(timeTreshold),' seconds'));
    xlabel('X axis');
    ylabel('Y axis');
    zlabel('Z axis');
    
end