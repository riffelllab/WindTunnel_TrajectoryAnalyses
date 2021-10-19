% Function to load the wind tunnel test section structure in the plot
% Returns a plot figure containing a representation of the test section
%volume
function figTS = load_test_section_volumen(cluesSetup)
    %test section X,Y,Z limits (start_x= -lim_x and end_x= +lim_x --> Total x distance is lim_x*2)
    lim_x= 0.9144;
    lim_y= 0.3048;
    %lim_z= 0.3048;
    lim_z= 0.6;
    fixed_x= linspace(-lim_x, lim_x, 3);
    fixed_y= linspace(-lim_y, lim_y, 3);
    fixed_z= linspace(0, lim_z, 3);
    
    auxX1(1:length(fixed_y))=fixed_x(1);
    auxX3(1:length(fixed_y))=fixed_x(3);
    auxY1(1:length(fixed_x))=fixed_y(1);
    auxY3(1:length(fixed_x))=fixed_y(3);
    auxZ1(1:length(fixed_x))=fixed_z(1);
    auxZ3(1:length(fixed_x))=fixed_z(3);

    test= [fixed_z(1),fixed_z(1),fixed_z(1)];
    figTS= figure;
   
    %load the borders of the test section
    plot3(auxX1, fixed_y, test, 'k'); 
    hold on
    plot3(auxX1, fixed_y, [fixed_z(3),fixed_z(3),fixed_z(3)], 'k');
    plot3(auxX3, fixed_y, [fixed_z(1),fixed_z(1),fixed_z(1)], 'k'); 
    plot3(auxX3, fixed_y, [fixed_z(3),fixed_z(3),fixed_z(3)], 'k');
    plot3(auxX1, [fixed_y(1),fixed_y(1),fixed_y(1)], fixed_z, 'k'); 
    plot3(auxX1, [fixed_y(3),fixed_y(3),fixed_y(3)], fixed_z, 'k');
    plot3(auxX3, [fixed_y(1),fixed_y(1),fixed_y(1)], fixed_z, 'k'); 
    plot3(auxX3, [fixed_y(3),fixed_y(3),fixed_y(3)], fixed_z, 'k');
    plot3(fixed_x, auxY1, auxZ1, 'k'); 
    plot3(fixed_x, auxY3, auxZ1, 'k'); 
    plot3(fixed_x, auxY1, auxZ3, 'k'); 
    plot3(fixed_x, auxY3, auxZ3, 'k'); 
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
    
    %Load the clues used (if any)
    if isempty(cluesSetup) == 0
        %Load the XY position for ODOR
        x= -lim_x + cell2mat(cluesSetup(1,2));
        y= cell2mat(cluesSetup(1,3));
        z= cell2mat(cluesSetup(1,4));
        plot3(x,y,z, 'p', 'markerfacecolor','magenta');
        for i= 2:length(cluesSetup(:,1))
            %Load the XY position for each visual clue
            x= -lim_x + cell2mat(cluesSetup(i,2));
            y= cell2mat(cluesSetup(i,3));
            plot(x, y, 'p', 'markerfacecolor', char(cluesSetup(i,1)));
        end;
    end;
    hold off
    %grid on