% clc; clear; close all;
year=2020;
doy=46;

load coastlines; % Load Earth coastline data
coastlat;
coastlon(coastlon<-120)=coastlon(coastlon<-120)+360; % Shift longitude for easier plotting
%% Plot 2D map (Sub-satellite points)
figure; hold on; grid on;
scatter(coastlon, coastlat,2, [1/3 1/3 1/3], 'filled'); % Plot the Earth map

%% Plot sub-satellite points of satellite trajectories
year=num2str(year);
% Paths
tracePath="C:\Users\90576\Desktop\FlexPower\Paper\plotCode\figMode\satData\"+year;% Trajectory files; each file contains the latitude and longitude of 32 satellites for one day
timePath="C:\Users\90576\Desktop\FlexPower\Paper\plotCode\figMode\abnormal_trace_time\"+year;% Anomaly time files
filename=timePath+"\1.csv";
% Read CSV file using readtable and automatically use headers as variable names
T = readtable(filename);  
% Get table headers (column names) 
columnNames = T.Properties.VariableNames;  
% Initialize an array to store processed column names
processedNames = zeros(1, length(columnNames));  
% Iterate through each column name, split it, and extract the part before "_"
for i = 1:length(columnNames)  
    % Split the column name 
    parts = strsplit(columnNames{i}, '_');  
    % Take the first part; if empty after splitting, keep the original column name  
    if ~isempty(parts)  
        processedNames(i) =str2num(parts{1}(2:4));  
    else  
        processedNames(i) = str2num(columnNames{i}(2:4));  
    end  
end  
doys=unique(processedNames); % The provided dataset may not include every day of the year
doyIndex = find(doys==doy);
str=[];
latEdge = []; lonEdge = [];
for PRN=1:32 % Number of satellites
    doy=doys(doyIndex);
    tracefilename=tracePath+"\"+num2str(doy)+".csv"; % File for the specific day
    try
        timefilename=timePath+"\"+num2str(PRN)+".result.csv";
        data=readmatrix(timefilename); % Every two columns represent one day
    catch
        timefilename=timePath+"\"+num2str(PRN)+".csv";
        data=readmatrix(timefilename);% Every two columns represent one day
    end    
    % Read CSV file
    index=2*doyIndex-1;
    if ~isnan(data(1,index))
        startRow=data(:,index);% Determine the start time of the abnormal trajectory
        startRow=startRow(~isnan(startRow));
        endRow=data(:,index+1);% Determine the end time of the abnormal trajectory
        endRow=endRow(~isnan(endRow));
        data = readmatrix(tracefilename);
        num=length(startRow);% Determine the number of trajectory segments    
        for i=1:num
            time=[startRow(i);endRow(i)];
            dataSubset = data(time, :);
            % Separate longitude and latitude data
            index=PRN*2-1;
            longitudes = dataSubset(:,index+1); % Longitude
            longitudes((longitudes<-120))=longitudes((longitudes<-120))+360; 
            latitudes = dataSubset(:,index);% Latitude
            lonEdge=[lonEdge;longitudes];
            latEdge=[latEdge;latitudes];   
        end
    end
    h=scatter(lonEdge,latEdge,20,'filled');
end

dpi = 1;
coneVote = {};
for ibatCone = 1:dpi:size(batCone,1)
    disp(ibatCone)
    for jbatCone =  1:dpi:size(batCone,2)
        coneVote{ibatCone,jbatCone} = [];
        for iEdge = 1:length(latEdge)
            latVote = abs(batCone{ibatCone,jbatCone}(1,:)-latEdge(iEdge))<2;
            lonVote = abs(batCone{ibatCone,jbatCone}(2,:)-lonEdge(iEdge))<2;
            if sum(latVote.*lonVote)>0
                coneVote{ibatCone,jbatCone} = [coneVote{ibatCone,jbatCone},[latEdge(iEdge),lonEdge(iEdge)]'];
            end
        end
    end
end

copyLatEdge = latEdge;
copyLonEdge = lonEdge;
centers = [];
maxSort = 1;
while length(copyLonEdge)>3
    lenConeVote = cellfun(@(x) size(x, 2), coneVote);
    [maxLen,maxLenIndex] = maxk(lenConeVote(:),maxSort);   
    [centersRow,centersCol] = ind2sub(size(lenConeVote),maxLenIndex(maxSort));
    centers = [centers,[centersRow,centersCol,maxLen(maxSort)]'];
    fixLatLon = coneVote{centersRow,centersCol};
    for ilenFix = 1:length(fixLatLon)
        disp(copyLonEdge')
        if isempty(copyLonEdge)
            break
        end
        delLatLon = fixLatLon(:,ilenFix);
        copyLatEdge(copyLatEdge==delLatLon(1)) = [];
        copyLonEdge(copyLonEdge==delLatLon(2)) = [];
        for iCone = 1:size(coneVote,1)
            for jCone =  1:size(coneVote,2)
                if (iCone == centersRow) && (jCone == centersCol)
                    continue                    
                end
                if ~isempty(coneVote{iCone,jCone})
                    coneVote{iCone,jCone}(coneVote{iCone,jCone}==delLatLon) = [];
                    if min(size(coneVote{iCone,jCone}))==1
                        coneVote{iCone,jCone} = reshape(coneVote{iCone,jCone},2,length(coneVote{iCone,jCone})/2);
                    end
                end
            end
        end
    end
    maxSort = maxSort + 1;
end



