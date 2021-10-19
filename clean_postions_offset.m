function data = clean_position_offset(data, lim_x, lim_y, lim_z)
    %find the positions outside the test section and delete them (they are
    %noise or wrong values, not to use)
    %data= Struct with the following fields:
    %   attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z
      
    IdxX = find(data.attr_x > lim_x | data.attr_x < -lim_x);
    IdxY = find(data.attr_y > lim_y | data.attr_y < -lim_y);
    if lim_z== 0.3048
        IdxZ = find(data.attr_z > lim_z | data.attr_z < -lim_z);
    else
        IdxZ = find(data.attr_z > lim_z | data.attr_z < 0);
    end;
    %IdxZ = find(posZ > 0.8 | posZ < -0.2);

    %Clean these values outside the test section
    WrongIdx = [IdxX; IdxY; IdxZ];
    data.attr_x(WrongIdx,:) = [];
    data.attr_y(WrongIdx,:) = [];
    data.attr_z(WrongIdx,:) = [];
end
