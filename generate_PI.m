
% Generate Preference Index (PI) between a base and color cues (values in col 1 and 2 from data)
% PI= (data(i)- data(k))/totalTraj
% Arguments:
%   - data: Array containing the values to calculate the PI (row= expIndex)
%   - baseCuesIndexList: list with the indexes assigned to the white
%                         visual cue (in dataset(:).expCues
% Output:
%   - piList: List with the PI towards the Main Value

function piList= generate_PI(data, baseCueIndexList)
    
    %data=countListWithCO2;
    piList= zeros(length(data(:,1)),1);
    for expIndex= 1:length(data(:,1))
        totalTraj= sum(data(expIndex,:));
        if baseCueIndexList(expIndex)== 1
            colorCueIndex=2;
        else
            colorCueIndex= 1;
        end
        pi= (data(expIndex,colorCueIndex)- data(expIndex,baseCueIndexList(expIndex)))/totalTraj;
        piList(expIndex)= pi;      
    end
    
    
end


% test
% duration= 120;
% p1BaseC= find(baseColorIndexList == 1);
% p2BaseC= find(baseColorIndexList== 2);
% %Do the same for the counts when the CO2 is being released
% sumCountsInBaseClrCO2= cumsum([p1CO2(p1BaseC, 1:duration); p2CO2(p2BaseC, 1:duration)], 2);
% sumCountsInTestClrCO2= cumsum([p2CO2(p1BaseC, 1:duration); p1CO2(p2BaseC, 1:duration)], 2);
% totalCountsCO2= sumCountsInBaseClrCO2(:,end) +sumCountsInTestClrCO2(:,end);
% piList3= bsxfun(@rdivide,(sumCountsInTestClrCO2(:,end)- sumCountsInBaseClrCO2(:,end)),totalCountsCO2);
% for i=1:4
%     pi= (sumCountsInTestClrCO2(i,end)- sumCountsInBaseClrCO2(i,end))/totalCountsCO2(i);
%     piList2(i,1)= pi;  
% end


% Generate Preference Index (PI) for values in col 1 and 2 from data 
% PI= (data(i)- data(k))/totalTraj
% Arguments:
%   - data: Array containing the values to calculate the PI (row= expIndex)
%   - i:    Index for the Main Value in the expCues list of first experiment
%   - k:    Index for the 2nd value in the expCues list of first experiment
% Output:
%   - piList: List with the PI towards the Main Value

% function piList= generate_PI(data, i, k)
%     
%     %data=countListWithCO2;
%     piList= zeros(length(data(:,1)),1);
%     for expIndex= 1:length(data(:,1))
%         whiteCueIndex= find(strcmp(cellstr(varTemp),'white'));
%         disp(strcat('blanquito: ',num2str(whiteCueIndex)));
%         totalTraj= sum(data(expIndex,:));
%         pi= (data(expIndex,i)- data(expIndex,k))/totalTraj;
%         piList(expIndex)= pi;
%         t=i;
%         i=k;
%         k=t;
%     end
% end