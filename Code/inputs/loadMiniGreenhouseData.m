function [weather, startTime] = loadMiniGreenhouseData(firstDay, seasonLength)
    SECONDS_IN_DAY = 24*60*60;

    currentFile = mfilename('fullpath');
    currentFolder = fileparts(currentFile);

    path = [currentFolder '\minigreenhouse_2.mat'];

    %% load hi res 
    minigreenhouse = load(path).datasetminigreenhouse;
    
    %% Cut out the required season
    interval = minigreenhouse(2,1) - minigreenhouse(1,1); % assumes all data is equally spaced

    firstDay = mod(firstDay, 365); % in case value is bigger than 365
    
    % use only needed dates
    startPoint = 1+round((firstDay-1)*SECONDS_IN_DAY/interval);
        % index in the time array where data should start reading
    endPoint = startPoint-1+round(seasonLength*SECONDS_IN_DAY/interval);

    % calculate date and time of first data point
    startTime = datetime(2000,1,1,0,0,0)+minigreenhouse(startPoint,1)/SECONDS_IN_DAY;
       
    dataLength = length(minigreenhouse(:,1));
    newYears = (endPoint-mod(endPoint,dataLength))/dataLength; 
        % number of times data crosses the new year
    
    if endPoint <= dataLength % required season passes over end of year
        season = minigreenhouse(startPoint:endPoint,:);
    else
        season = minigreenhouse(startPoint:end,:);
        for n=1:newYears-1
            season = [season; minigreenhouse];
        end
        endPoint = mod(endPoint, dataLength);
        season = [season; minigreenhouse(1:endPoint,:)];
    end
    % Reformat data
    %% INCREASE TEMPERATURE BY 1.5°C, TO FIT BETTER WITH MODERN CLIMATE!!!
    %   length  - length of desired dataset (days)
    %   weather  will be a matrix with 9 columns, in the following format:
    %       weather(:,1)    timestamps of the input [datenum] in 5 minute intervals
    %       weather(:,2)    radiation     [W m^{-2}]  outdoor global irradiation 
    %       weather(:,3)    temperature   [°C]        outdoor air temperature
    %       weather(:,4)    humidity      [kg m^{-3}] outdoor vapor concentration
    %       weather(:,5)    co2 [kg{CO2} m^{-3}{air}] outdoor CO2 concentration
    %       weather(:,6)    wind        [m s^{-1}] outdoor wind speed
    %       weather(:,7)    sky temperature [°C]
    %       weather(:,8)    temperature of external soil layer [°C]
    %       weather(:,9)    daily radiation sum [MJ m^{-2} day^{-1}]

    weather(:,1) = interval*(0:length(season(:,1))-1); % time
    weather(:,2) = season(:,5); % lux
    weather(:,3) = season(:,3) + 1.5; % temperature
    weather(:,4) = rh2vaporDens(weather(:,3), season(:,4)); % humidity not in PPM
    weather(:,5) = co2ppm2dens(weather(:,3), season(:,2));
    % weather(:,6) = 1*ones(length*288,1);
    weather(:,6) = 0;
    weather(:,7) = weather(:,3) - 20;
    % weather(:,8) = 20*ones(length*288,1);
    
    % convert timestamps to datenum
    % weather(:,1) = time/86400+1;
    % weather(:,9) = dayLightSum(weather(:,1), weather(:,2));

end


