% Function that returns the MATLAB color code assigned to a given color name. 
% - Arguments
%       - color: Name or code of the color used as visual cue in the experiment
% - Returns
%       - colorRGB: RGB representation of the color used as visual cue.
%           This value will used when plotting the experiment setup
function colorRGB = convert_color(color)
    switch(char(color))
        case 'gray9.5' 
            % Almost White (old gray1)
            colorRGB= [0, 0, 0]+0.95;
        case 'gray6.5' 
            % (old gray2)
            colorRGB= [0, 0, 0]+0.65;
        case 'gray4.5' 
            % (old gray3)
            colorRGB= [0, 0, 0]+0.45;
        case 'gray2.5'
            % Almost black (old gray4)
            colorRGB= [0, 0, 0]+0.25;
        case 'gray4.0'
            % Gray used to compare with the R-HUE red 
            %colorRGB= [0, 0, 0]+0.18;
            colorRGB= [0, 0, 0]+0.40;
        case 'RW-HUE'
            % ColorAid first Red
            colorRGB= [0.8500, 0.3250, 0.0980];
        case 'R-HUE'
            % ColorAid second Red
            colorRGB=[0.6350, 0.0780, 0.1840];
        case 'RO-S1'
            % ColorAid third Red
            colorRGB= [1 0 1];
        case 'R-T1'
            % ColorAid fourth Red
            colorRGB= [0.40 0 0.20];
        case 'RS1'
            colorRGB= [0.60 0 0.20];
        case 'RC-HUE'
            colorRGB= [0.80 0 0.20];
        case 'YwHue'
            % ColorAid Yellow
            colorRGB= [1, 1, 0];
        case 'BwT1'
            % ColorAid Blue
            colorRGB= [0, 0, 1];
        case 'BVT2'
            % ColorAid Purple
            colorRGB= [0.5, 0, 0.5];
        case 'GwT1'
            % ColorAid first Green
            colorRGB= [0, 1, 0];
        case 'YGcHue'
            % ColorAid second Green 
            colorRGB= [0, 0.6, 0];
        case 'GcT1'
            % ColorAid third Green
            colorRGB= [0, 0.3, 0];
        case 'GS1'
            % ColorAid 4th Green
            colorRGB= [0, 0.1, 0];
        case 'GwT2'
            % ColorAid 5th Green
            colorRGB= [0.2, 1, 0.2];
        otherwise
            colorRGB= char(color);
    end
end