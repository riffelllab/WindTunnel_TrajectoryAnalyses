
%Function to load the insect IDs, their XYZ values, frames number and timestamps when the
%position was acquired
%Arguments: 
%   - File: Path/name of the h5 file
%   - LoadFullDataset: Boolean flag to check if we want to work with the
%   - Flag to select full dataset or only the first half
%Returns: object ids, frame numbers and X,Y,Z positions

function [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z]= load_data_from_file(file, loadFullDataset)

    %location_ML='/ML_estimates';   %dataset name where the attributes are contained
    location_KA='/kalman_estimates';   %dataset name where the attributes are contained

    % from the DATASET ML_estimates and kalman_estimates 
    %I will read the attributes frame, x y and z to plot them
    dataset= h5read(file,location_KA); 

    if loadFullDataset == true
        %Load the full dataset
        disp(strcat(' * Loading full dataset from ',{' '},file(end-27:end-3)));
        attr_frame=double(dataset.frame);
        attr_time= double(dataset.timestamp); 
        attr_x=double(dataset.x);
        attr_y=double(dataset.y);
        attr_z=double(dataset.z);     
        %attr_id= dataset.obj_id(:);
        attr_id= double(dataset.obj_id);
    else
        %Load only the first half of the dataset
        disp(strcat(' * Loading only the first half of the dataset from ',{' '},file(end-27:end-1)));
        halfSize= ceil(length(dataset.frame)/2);   %load half dataset
        attr_frame=dataset.frame(1:halfSize);
        attr_time= dataset.timestamp(1:halfSize); 
        attr_x=dataset.x(1:halfSize);
        attr_y=dataset.y(1:halfSize);
        attr_z=dataset.z(1:halfSize);     
        attr_id= dataset.obj_id(1:halfSize);
    end
end