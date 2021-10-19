% Generates the output folder in function of the input folder name and
% return the name of the output folder
% Arguments:
%   - outputPath: path where to create the experiment output folder
%   - inputFolder: Name of the experiment input folder
% Output
%   - outFolder: The name of the output folder generated

function outputFolder = load_output_folder(outputPath, inputFolder)

    %Header used in inputFolder and outputFolder
    inputHeader= 'FLYDRA_Trajectories_';
    outputHeader= 'exp_';
    outputDataFolder= 'dataFolder\';
    %Name of the output folder
    outputFolder= strcat(outputHeader, erase(inputFolder, inputHeader));
    
    [status, msg, msgId] = mkdir(outputPath, outputFolder);
    if ~status
        disp('  - WARNING: Output folder NOT created');
            % Create in the output folder, the folder for the .xls/.csv files
            [status, msg, msgId] = mkdir(strcat(outputPath, outputFolder), outputDataFolder );
    else
        disp(newline);
        disp(strcat(' * Folder ',{' '},outputFolder, ' created or already exist.'));
        disp(strcat(' * Folder path= ',{' '},outputPath, newline));
        disp(' ====================   ==================== ');
    end
end