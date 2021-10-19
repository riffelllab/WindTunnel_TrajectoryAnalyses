
%function to plot the deltaX, deltaY and deltaY between frames
%(how big is the step between frames for each axis)
%Arguments: 
%   - attr_x: List containing all the positions in the X axis tracked by Flydra
%   - attr_y: List containing all the positions in the Y axis tracked by Flydra
%   - attr_Z: List containing all the positions in the X axis tracked by Flydra
%Returns: Nothing. a Plot is created showing the delta values

%WARNING: - Some variations values are bigger than the limits set in the Y axis axis, 
%         but the yLim() allows us to see the small changes. 
%         - To see the full variation between positions, comment the ylim()
%         lines

function plot_step_XYZ_sizes(attr_x, attr_y, attr_z)
    axisLimit= 1.5;
    diffX= diff(attr_x);
    diffY= diff(attr_y);
    diffZ= diff(attr_z);
    figure(3)
    subplot(3,1,1);
    plot(diffX,'b');
    title('Plot step variation for X axis');
    xlabel('Frame number');
    ylim([-axisLimit axisLimit]);
    ylabel('step variation');
    subplot(3,1,2);
    plot(diffY,'r');
    title('Plot step variation for Y axis');
    xlabel('Frame number');
    ylim([-axisLimit axisLimit]);
    ylabel('step variation');
    subplot(3,1,3);
    plot(diffZ,'g');
    title('Plot step variation for Z axis');
    xlabel('Frame number');
    ylim([-axisLimit axisLimit]);
    
    ylabel('step variation');
    disp(' * Plot Y axis have been limited between -5 and 5.')
    disp(strcat('   - Biggest step variation between frames in the X axis: ',num2str(max(diffX))))
    disp(strcat('   - Biggest step variation between frames in the Y axis: ',num2str(max(diffY))))
    disp(strcat('   - Biggest step variation between frames in the Z axis: ',num2str(max(diffZ))))

    disp(strcat('   - Mean step variation between frames in the X axis: ',num2str(mean(diffX))))
    disp(strcat('   - Mean step variation between frames in the Y axis: ',num2str(mean(diffY))))
    disp(strcat('   - Mean step variation between frames in the Z axis: ',num2str(mean(diffZ))))
end
