
% Function to plot the heatmaps containing the XYZ preferences.
% The data is displayed in 2 subplots: XY positions and XZ positions.
% -ve Y axis is the left side of the WT and +ve Y axis is the right
% If required, odor and Visual clues are also plotted
% The graph is saved in the OUTPUT folder assigned for the experiment
% (information generated in the base workspace, not in this fct)
% Arguments:
%   - plotCues: Control falg to add the odor and visual cues used to the plot
%   - posX: column vector with the positions in the X axis 
%   - posY: column vector with the positions in the Yaxis 
%   - posZ: column vector with the positions in the Z axis
%   - cuesSetup: cell matrix containing the cue and its XYZ position 
%               [{cue} {posX} {posY} {posZ}]
%   - imgTitle: name to use when saving the heatmap as .png or .fig
function plot_XY_XZ_heatmaps_v5(plotCues, posX, posY, posZ, cuesSetup, topValueNormalized, imgTitle)
    disp(imgTitle)
    switch (imgTitle(end-6:end))
        case 'prevCO2'
            t1= ('Heatmap intial air only');
        case 'withCO2'
            t1= 'Heatmap air + CO2';
        case 'postCO2'
            t1= 'Heatmap final air only';
        otherwise
            t1= 'Full_dataset';
    end
    
    %topValueNormalized= 0.0001; %0.0003;
    nBins= [200,600];
    
    
    t2= strcat(' - Data normalized (max=', num2str(topValueNormalized),') -');
    t3= '(Top view)';
    titleText= sprintf(strcat(t1,'\n',t2,'\n',t3));

    figure(gcf)
    subplot(2,1,1);
    %h= histogram2(posX, posY,[600 200],'DisplayStyle','tile','ShowEmptyBins','on');
    %h= histogram2(posX, posY, [200 600],'Normalization', 'probability', 'DisplayStyle','tile','ShowEmptyBins','on');
    h= histogram2(posX, posY, nBins,'Normalization', 'probability', 'DisplayStyle','tile','ShowEmptyBins','on');
    %valuesUpdated = h.Values/sum(h.Values(:));
    colormap(jet);
    %imagesc(h.Values');
    %imagesc(valuesUpdated*100);
    colorbar();
    caxis([0 topValueNormalized]);
    
    %set(gca, 'YDir','reverse')  %To keep the Y axis as FLydra  consider it (-ve values close towards the white sheet side in the wind tunnel
    xlabel('X')
    ylabel('Y')
    title(titleText);
    
    % Add the odor and visual cues to the heatmap plot (color and XZ position)
    %h= add_visual_cues_to_plot(h, cuesSetup(:,[1 2 3]), plotCues);

    
    subplot(2,1,2);
    %h= histogram2(posX, posZ, [600 200], 'DisplayStyle','tile','ShowEmptyBins','on');
    h2= histogram2(posX, posZ, nBins, 'Normalization', 'probability', 'DisplayStyle','tile','ShowEmptyBins','on');
    %valuesUpdated = h.BinCounts/sum(h.BinCounts(:));
    colormap(jet);
    %imagesc(h.Values');
    %imagesc(valuesUpdated*100);
    colorbar();
    caxis([0 topValueNormalized]);
    xlabel('X')
    ylabel('Z')
    title('(Side view)')
    
    % if plotCues= true, Add the odor and visual cues to the heatmap plot
    % (color and XZ position)
    if plotCues
        subplot(2,1,1);
        h= add_visual_cues_to_plot(h, cuesSetup(:,[1 2 3]));
        subplot(2,1,2);
        h2= add_visual_cues_to_plot(h2, cuesSetup(:,[1 2 4]));
        imgTitle= strcat(imgTitle, '_vc');
    else
        %If we aren't plotting the cues, hide the axis too 
        subplot(2,1,1);
        set(gca,'visible','off')
        subplot(2,1,2);
        set(gca,'visible','off')
        %set(gcf, 'xticklabel', []);
        %set(gcf,'XColor', 'none','YColor','none')
    end
    
    % Save image
    save_plot_in_exp_folder(gcf, imgTitle);
%     outputPath = evalin('base', 'outputPath');
%     outputFolder = evalin('base', 'outputFolder'); 
%     disp('path teasted')
%     disp(strcat(outputPath, outputFolder))
    %print(gcf,strcat(outputPath, outputFolder, imgTitle, '.png'),'-dpng','-r300');        % *// 300 dpi
    %print(gcf,strcat(imgTitle, '.tiff'),'-dtiff','-r300');        % *// 300 dpi
end