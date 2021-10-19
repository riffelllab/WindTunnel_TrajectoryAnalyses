% Function to save the conts inside a volume over time (counts occurred at
% a given number of secods)
% Arguments:
%   - pX: Matrix with the number of counts in the volume located in position X 
%         for each time group per experiment. (rows: experiment index, cols: time group)
function save_counts_over_time_in_file(p1, p2, p3, p4)

    % add position value to the matrices
    p1= horzcat(ones(length(p1(:,1)), 1),p1);
    p2= horzcat(ones(length(p2(:,1)),1)*2,p2);
    % group the values for position 1 and 2
    data= vertcat(p1, p2);

    %If we are using experiments with 4 cues
    if (nnz(p3) && nnz(p4))
        % add position value to the matrices
        p3= horzcat(ones(length(p3(:,1)), 1)*3,p1);
        p4= horzcat(ones(length(p4(:,1)),1)*4,p2);
        % Add the values for position 3 and 4 to the group with pos 1 and 2
        data= vertcat(data, p3);
        data= vertcat(data, p4);
    end

    outputPath = evalin('base', 'outputPath');
    outputFolder = evalin('base', 'outputFolder');  
    path=strcat(outputPath, outputFolder);
    % Columns in the csv file: Position, time grp1, time grp2, ..., time grp N
    % Each row represent 1 experiment
    writematrix(data,strcat(path, 'countsPerTimeGrp.csv')); 

end