function [outdoor, indoor, controls, startTime] = loadMiniGreenhouseData2(firstDay, seasonLength)
% loadMiniGreenhouseData2 Get data from real mini-greenhouse experiments
% The following datasets are available:
% - 
% - 
% The data is given in 5-minute intervals.
% Based on:
% David Katzin, Simon van Mourik, Frank Kempkes, and Eldert J. Van Henten. 2020. 
% "GreenLight - An Open Source Model for Greenhouses with Supplemental Lighting: Evaluation of Heat Requirements under LED and HPS Lamps.” 
% Biosystems Engineering 194: 61–81. https://doi.org/10.1016/j.biosystemseng.2020.03.010
% 
% 
% Usage:
%   [outdoor, indoor, contorls, startTime] = loadGreenhouseData(firstDay, seasonLength, type)
% The dataset contain a table in the following format:
% Column    Description                         Unit             
% 1 		Time 								datenum 
% 2 		Radiation global				    W m^{-2} outdoor global irradiation 
% 3         Radiation inside                    W m^{-2}
% 4 		Temp in 							°C
% 5         Temp out                            °C
% 6 		Relative humidity in 				%	
% 7         Relatve humidity out                %
% 8         CO2 in                              ppm
% 9         CO2 out                             ppm
% 10        Toplights on/off                    0/1 (1 is on)
% 11        Average roof ventilation aperture	(average between lee side and wind side)	0-1 (1 is fully open)
%
% Output:
% Function inputs:
%   lampType        Type of lamps in the greenhouse. Choose between 
%                   'hps', 'led', or 'none' (default is none)
%   weather         A matrix with 8 columns, in the following format:
%       weather(:,1)    timestamps of the input [s] in regular intervals
%       weather(:,2)    radiation     [W m^{-2}]  outdoor global irradiation 
%       weather(:,3)    temperature   [°C]        outdoor air temperature
%       weather(:,4)    humidity      [kg m^{-3}] outdoor vapor concentration
%       weather(:,5)    co2 [kg{CO2} m^{-3}{air}] outdoor CO2 concentration
%       weather(:,6)    wind        [m s^{-1}] outdoor wind speed
%       weather(:,7)    sky temperature [°C]
%       weather(:,8)    temperature of external soil layer [°C]
%       weather(:,9)    daily radiation sum [MJ m^{-2} day^{-1}]
%   startTime       date and time of starting point (datetime)
%
%   controls        (optional) A matrix with 8 columns, in the following format:
%       controls(:,1)     timestamps of the input [s] in regular intervals of 300, starting with 0
%       controls(:,2)     Energy screen closure 			0-1 (1 is fully closed)
%       controls(:,3)     Black out screen closure			0-1 (1 is fully closed)
%       controls(:,4)     Average roof ventilation aperture	(average between lee side and wind side)	0-1 (1 is fully open)
%       controls(:,5)     Pipe rail temperature 			°C
%       controls(:,6)     Grow pipes temperature 			°C
%       controls(:,7)     Toplights on/off                  0/1 (1 is on)
%       controls(:,8)     Interlight on/off                 0/1 (1 is on)
%       controls(:,9)     CO2 injection                     0/1 (1 is on)
%
%   indoor          (optional) A 3 column matrix with:
%       indoor(:,1)     timestamps of the input [s] in regular intervals of 300, starting with 0
%       indoor(:,2)     temperature       [°C]             indoor air temperature
%       indoor(:,3)     vapor pressure    [Pa]             indoor vapor concentration
%       indoor(:,4)     co2 concentration [mg m^{-3}]      indoor co2 concentration

    SECONDS_IN_DAY = 24*60*60;
    
    %% load file
    currentFile = mfilename('fullpath');
    currentFolder = fileparts(currentFile);
    
    path = [currentFolder '\minigreenhouse_3.mat'];
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
    % startTime = datetime(2000,1,1,0,0,0)+minigreenhouse(startPoint,1)/SECONDS_IN_DAY;
    startTime = datetime(minigreenhouse(1,1),'ConvertFrom','datenum');

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
    %% REFORMAT DATA
    % See the format of the dataset above

    %% WEATHER
    % Weather reformartted dataset
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

    outdoor(:,1) = interval*(0:length(season(:,1))-1); % time
    outdoor(:,2) = season(:,2); 
    outdoor(:,3) = season(:,5);
    outdoor(:,4) = rh2vaporDens(outdoor(:,3), season(:,7)); % Convert relative humidity [%] to vapor density [kg{H2O} m^{-3}]
    outdoor(:,5) = co2ppm2dens(outdoor(:,3), season(:,7)); % Convert CO2 molar concetration [ppm] to density [kg m^{-3}]
    outdoor(:,6) = 0;
    

    %% INDOOR
    % Indoor reformartted dataset
    %   indoor          (optional) A 3 column matrix with:
    %   indoor(:,1)     timestamps of the input [s] in regular intervals of 300, starting with 0
    %   indoor(:,2)     temperature       [°C]             indoor air temperature
    %   indoor(:,3)     vapor pressure    [Pa]             indoor vapor concentration
    %   indoor(:,4)     co2 concentration [mg m^{-3}]      indoor co2 concentration
    indoor(:,1) = outdoor(:,1);
    indoor(:,2) = season(:,4);
    indoor(:,3) = 0;
    indoor(:,4) = 1e6*co2ppm2dens(indoor(:,2), season(:,9)); % convert co2 from ppm to mg m^{-3}
    
    %% CONTROL
    % Control reformartted dataset
    %   controls        (optional) A matrix with 8 columns, in the following format:
    %   controls(:,1)     timestamps of the input [s] in regular intervals of 300, starting with 0
    %   controls(:,2)     Energy screen closure 			0-1 (1 is fully closed)
    %   controls(:,3)     Black out screen closure			0-1 (1 is fully closed)
    %   controls(:,4)     Average roof ventilation aperture	(average between lee side and wind side)	0-1 (1 is fully open)
    %   controls(:,5)     Pipe rail temperature 			°C
    %   controls(:,6)     Grow pipes temperature 			°C
    %   controls(:,7)     Toplights on/off                  0/1 (1 is on)
    %   controls(:,8)     Interlight on/off                 0/1 (1 is on)
    %   controls(:,9)     CO2 injection                     0/1 (1 is on)
    controls(:,1) = outdoor(:,1);
    controls(:,2) = 0;
    controls(:,3) = 0; 
    controls(:,4) = season(:,11);
    controls(:,5) = indoor(:,2);        % pipe rail temperature with adjustment for air
    controls(:,6) = indoor(:,2);        % grow pipes temperature with adjustment for air
    controls(:,7) = indoor(:,10);
    controls(:,8) = 0;
    controls(:,9) = 0;

end