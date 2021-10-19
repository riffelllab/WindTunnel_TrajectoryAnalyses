%Function to control if the path tracked is inside the test section or not.
%If trajectory is outside test section it can be due to an artifact (dust,
%insect outside the Wind Tunnel,...) or that the calibration is not valid
%anymore.
%Arguments: 
%   - objID: Insect ID associated to the current insect trajectory
%   - objXYZ: Current insect trajectory
%   - lim_x, lim_y and lim_z: current limits for the X,Y and Z axis. These
%   limits are the limits of the test section volume


function controlFlag= check_traj_inside(objID, objXYZ, lim_x, lim_y, lim_z)
    lim_x= 0.9144;
    lim_y= 0.3048;
    lim_z= 0.6096;
    Ix= find(objXYZ(:,1) < lim_x & objXYZ(:,2) > -lim_x);
    Iy= find(objXYZ(:,2) < lim_y & objXYZ(:,2) > -lim_y);         
    %Iz= find(objXYZ(:,3) < lim_z & objXYZ(:,3) > -lim_z);
    Iz= find(objXYZ(:,3) < lim_z & objXYZ(:,3) > 0);
    if Ix > 0 | Iy > 0 | Iz> 0
        %disp(strcat(' * Outside the Test Section, ObjID: ',num2str(objID)));
        controlFlag= false;
    else
        disp(strcat(' * Inside the Test Section, ObjID: ',num2str(objID)));
        controlFlag= true;
    end;