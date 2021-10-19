% Function to erase the points detected outside the test section of the windtunnel (dust, error in the XYZ acquisitions, etc)
% TO USE IF WORKING WITH POSITIONS ONLY! 
% If working with POSITIONS AND SPEED, use erase_pos_outside_WT_v3
function  [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z] = erase_pos_outside_WT_v2(attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z, lim_x, lim_y, lim_z)
    %find the positions outside the test section and delete them (they are
    %noise or wrong values, not to use)
    data= [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z];

    IdxX = find(data(:,4) > lim_x | data(:,4) < -lim_x);
    IdxY = find(data(:,5) > lim_y | data(:,5) < -lim_y);
    if lim_z== 0.3048
        IdxZ = find(data(:,6) > lim_z | data(:,6) < -lim_z);
    else
        IdxZ = find(data(:,6) > lim_z | data(:,6) < 0);
    end;
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
end



% function [air, co2, postCo2] = erase_pos_outside_WT(air, co2, postCo2, lim_x, lim_y, lim_z)
%     %find the positions outside the test section and delete them (they are
%     %noise or wrong values, not to use)
%     %data= [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z];
%     
%     IdxX = find(air(:,4) > lim_x | air(:,4) < -lim_x);
%     IdxY = find(air(:,5) > lim_y | air(:,5) < -lim_y);
%     if lim_z== 0.3048
%         IdxZ = find(air(:,6) > lim_z | air(:,6) < -lim_z);
%     else
%         IdxZ = find(air(:,6) > lim_z | air(:,6) < 0);
%     end;
%     WrongIdx = [IdxX; IdxY; IdxZ];
%     
%     
%     IdxX = [find(air(:,4) > lim_x | air(:,4) < -lim_x); find(co2(:,4) > lim_x | co2(:,4) < -lim_x); find(postCo2(:,4) > lim_x | postCo2(:,4) < -lim_x)];
%     IdxY = [find(air(:,5) > lim_y | air(:,5) < -lim_y); find(co2(:,5) > lim_y | co2(:,5) < -lim_y); find(postCo2(:,5) > lim_y | postCo2(:,4) < -lim_y)];
%     if lim_z== 0.3048
%         IdxZ = [find(air(:,6) > lim_z | air(:,6) < -lim_z); find(co2(:,6) > lim_z | co2(:,6) < -lim_z); find(postCo2(:,6) > lim_z | postCo2(:,6) < -lim_z)];
%     else
%         IdxZ = [find(air(:,6) > lim_z | air(:,6) < 0); find(co2(:,6) > lim_z | co2(:,6) < 0); find(postCo2(:,6) > lim_z | postCo2(:,6) < 0)];
%     end;
%     %IdxZ = find(posZ > 0.8 | posZ < -0.2);
% 
%     %Clean these values outside the test section
%     WrongIdx = [IdxX; IdxY; IdxZ];
%     for value = WrongIdx
%         if value <= length(air)
%             air(value)=[];
%         elseif value > le
%     end;
%         
% end



% function [posX, posY, posZ] =erase_pos_outside_WT(posX, posY, posZ, lim_x, lim_y, lim_z)
%     %find the positions outside the test section and delete them (they are
%     %noise or wrong values, not to use)   
%     IdxX = find(posX > lim_x | posX < -lim_x);
%     IdxY = find(posY > lim_y | posY < -lim_y);
%     if lim_z== 0.3048
%         IdxZ = find(posZ > lim_z | posZ < -lim_z);
%     else
%         IdxZ = find(posZ > lim_z | posZ < 0);
%     end;
%     %IdxZ = find(posZ > 0.8 | posZ < -0.2);
% 
%     %Clean these values outside the test section
%     WrongIdx = [IdxX; IdxY; IdxZ];
%     posX(WrongIdx) = [];
%     posY(WrongIdx) = [];
%     posZ(WrongIdx) = [];
% end