% Function to save the plot in the output folder for the experiment type
% Arguments: 
%   - graph: plot to save as png
%   - imgTitle: name for the png file to create
function save_plot_in_exp_folder(graph, imgTitle)
    %Check if the output path and folder have been specified in the
    %workspace
    if (~isempty(evalin('base', 'outputPath')) &&  ~isempty(evalin('base', 'outputFolder')))
        outputPath = evalin('base', 'outputPath');
        outputFolder = evalin('base', 'outputFolder');  
        print(graph,strcat(outputPath, outputFolder, imgTitle, '.png'),'-dpng','-r300');        % *// 300 dpi
    else
        disp('  * Warning: output path and/or folder not specified in the variabes workspace ');
        disp('  * Plot saved in the AnalysisFLYDRA directory');
        print(graph, strcat(imgTitle, '.png'),'-dpng','-r300');        % *// 300 dpi
    
    end
end