% Function to erase the points detected outside the test section of the windtunnel (dust, error in the XYZ acquisitions, etc)
% TO USE IF WORKING WITH POSITIONS AND SPEED! 
% If working with POSITIONS ONLY, use erase_pos_outside_WT_v3
function  [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z, attr_xvel, attr_yvel, attr_zvel] = erase_pos_outside_WT_v3(attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z, attr_xvel, attr_yvel, attr_zvel, lim_x, lim_y, lim_z)
    %find the positions outside the test section and delete them (they are
    %noise or wrong values, not to use)
    data= [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z, attr_xvel, attr_yvel, attr_zvel];

    IdxX = find(data(:,4) > lim_x | data(:,4) < -lim_x);
    IdxY = find(data(:,5) > lim_y | data(:,5) < -lim_y);
    if lim_z== 0.3048
        IdxZ = find(data(:,6) > lim_z | data(:,6) < -lim_z);
    else
        IdxZ = find(data(:,6) > lim_z | data(:,6) < 0);
    end
    %IdxZ = find(posZ > 0.8 | posZ < -0.2);

    %Clean these values outside the test section
    WrongIdx = [IdxX; IdxY; IdxZ];
    data(WrongIdx,:) = [];
    
    %Reassign the values to the returned data
    attr_id=    data(:,1);
    attr_time=  data(:,2);
    attr_frame= data(:,3);
    attr_x=     data(:,4);
    attr_y=     data(:,5);
    attr_z=     data(:,6);
    attr_xvel=  data(:,7);
    attr_yvel=  data(:,8);
    attr_zvel=  data(:,9);
end


