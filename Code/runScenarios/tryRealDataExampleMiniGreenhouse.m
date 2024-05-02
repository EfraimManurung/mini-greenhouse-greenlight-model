% exampleMiniGreenhouse to run the GreenLight simulation
% Using createGreenLightModel

% Efraim Manurung, Information Technology Group
% Wageningen University
% efraim.efraimpartoginahotasi@wur.nl
% efraim.manurung@gmail.com

tic; % start the timer

% Weather argument for createGreenLightModel
seasonLength = 2; % season length in days
firstDay = 4; % days since beginning of data 

%% Choice of lamp
lampType = 'led'; % 'led', 'hps', or 'none'

% Interface
use_controls = true;
use_artificial_weather_dataset = false;
use_artificial_indoor_dataset = false;

% Controls dataset for createGreenLightModel instance
if use_controls
    controls = controlsParameters(seasonLength);
    controls(:,1) = (controls(:,1)-controls(1,1))*86400;
end

% Weather dataset for createGreenLightModel instance
if use_artificial_weather_dataset
    weather = weatherDataset(seasonLength);

    % convert weather timestamps from datenum to seconds since beginning of data
    weather(:,1) = (weather(:,1)-weather(1,1))*86400;
else
    % Real dataset for weather
    [weather, startTime] = loadMiniGreenhouseData(firstDay, seasonLength);
    % [weather, startTime] = loadSelYearHiRes(firstDay, seasonLength);

    secsInYear = seconds(startTime-datetime(year(startTime),1,1,0,0,0));
    % number of seconds since beginning of year to startTime

    weather(:,8) = soilTempNl(secsInYear+weather(:,1)); % add soil temperature
end

% Indoor dataset for createGreenLightModel instance
if use_artificial_indoor_dataset
    indoor = indoorDataset(seasonLength);
    indoor(:,1) = (indoor(:,1)-indoor(1,1))*86400;
end

%% Create an instance of createGreenLight with the default Vanthoor parameters
% conver weather timestamps from datenum to seconds since beginning of data
startTime = datetime(weather(1,1),'ConvertFrom','datenum');

% led = createGreenLightModel('led', weather, startTime, controls, indoor);
led = createGreenLightModel(lampType, weather, startTime);
led_controls = createGreenLightModel(lampType, weather, startTime,controls);

% Parameters for mini-greenhouse
setParamsMiniGreenhouse(led);      % set greenhouse structure
setMiniGreenhouseLedParams(led);   % set lamp params

setParamsMiniGreenhouse(led_controls);      % set greenhouse structure
setMiniGreenhouseLedParams(led_controls);   % set lamp params
% setParams4haWorldComparison(led); % set parameters according to a modern greenhouse
% setParams4haWorldComparison(led_controls); % set parameters according to a modern greenhouse

% Parameters for Bleiswjik2010, to make some tests
% setParamsBleiswijk2010(led);    % set greenhouse structure
% setBleiswijk2010LedParams(led); % set lamp params

%% Control parameters
setParam(led, 'rhMax', 87);        % upper bound on relative humidity  
setParam(led, 'lampsOn', 12);                     % Time of day when lamps go on                                                                  [hour]                                  0
setParam(led, 'lampsOff', 12);                   % Time of day when lamps go off
%setParam(led_controls, 'rhMax', 30);        % upper bound on relative humidity   

% setParam(led, 'thetaLampMax', 33.33);   %OK  % Maximum intensity of lamps

% Reset other parameters that depend on previously defined parameters
% setDepParams(led);
setDepParams(led_controls);

% Set initial values for crop
% led.x.cLeaf.val = 0.7*6240*10;
% led.x.cStem.val = 0.25*6240*10;
% led.x.cFruit.val = 0.05*6240*10;
% 
% % Set initial values for crop
% led_controls.x.cLeaf.val = 0.7*6240*10;
% led_controls.x.cStem.val = 0.25*6240*10;
% led_controls.x.cFruit.val = 0.05*6240*10;

%% Run simulation
solveFromFile(led, 'ode15s');
solveFromFile(led_controls, 'ode15s');

% set data to a fixed step size (5 minutes)
led = changeRes(led, 300);
led_controls = changeRes(led_controls, 300);

toc;
%% Plot some outputs 
% see setGlAux, setGlStates, setGlInput to see more options

% dateFormat = 'HH:00'; 
% This format can be changed, see help file for MATLAB function datestr

%% Figure 1 TEMPERATURE
figure(1)
plot(led.x.tAir,'LineWidth',1.5)
hold on
plot(led_controls.x.tAir,'LineWidth',1.5)
plot(led.d.tOut,'LineWidth',1.5)
hold off

% Get current x-axis ticks
numticks = get(gca,'XTick');

% Convert timestamps to datenum format
% divided by 86400
dateticks = datenum(datenum(led.t.label) + numticks / (60*60*24)); % Assuming timestamps are in seconds

% Format the date strings
if seasonLength > 2
    datestrings = datestr(dateticks, 'dd');
else
    datestrings = datestr(dateticks, 'HH:00');
end

% Set x-axis tick labels
xticks(numticks);
xticklabels(datestrings);

% Rotate x-axis tick labels to avoid overlapping
xtickangle(45);

% Add other plot elements
ylabel('Temperature (°C)')
%xlabel('Date')
xlabel('Time')
legend('Indoor','Indoor-controls','Outdoor')

%% Figure 2 RELATIVE HUMIDITY
figure(2)
plot(led.a.rhIn,'LineWidth',1.5)
hold on
plot(led_controls.a.rhIn, 'LineWidth',1.5)
plot(100*vp2dens(led.d.tOut,led.d.vpOut)./rh2vaporDens(led.d.tOut,100),'LineWidth',1.5);
hold off

% Get current x-axis ticks
numticks = get(gca,'XTick');

% Convert timestamps to datenum format
% divided by 86400
dateticks = datenum(datenum(led.t.label) + numticks / (60*60*24)); % Assuming timestamps are in seconds

% Format the date strings
if seasonLength > 2
    datestrings = datestr(dateticks, 'dd');
else
    datestrings = datestr(dateticks, 'HH:00');
end

% Set x-axis tick labels
xticks(numticks);
xticklabels(datestrings);

% Rotate x-axis tick labels to avoid overlapping
xtickangle(45);

% Add other plot elements
% xlabel('Date')
xlabel ('Time')
ylabel('Relative humidity (%)')
legend('Indoor', 'Indoor-controls','Outdoor')

%% Figure 3 CO2 IN PPM
figure(3)
plot(led.a.co2InPpm,'LineWidth',1.5)
hold on
plot(led_controls.a.co2InPpm,'LineWidth',1.5)
plot(co2dens2ppm(led.d.tOut,1e-6*led.d.co2Out),'LineWidth',1.5)
hold off

% Get current x-axis ticks
numticks = get(gca,'XTick');

% Convert timestamps to datenum format
% divided by 86400
dateticks = datenum(datenum(led.t.label) + numticks / (60*60*24)); % Assuming timestamps are in seconds

% Format the date strings
if seasonLength > 2
    datestrings = datestr(dateticks, 'dd');
else
    datestrings = datestr(dateticks, 'HH:00');
end

% Set x-axis tick labels
xticks(numticks);
xticklabels(datestrings);

% Rotate x-axis tick labels to avoid overlapping
xtickangle(45);

% Add other plot elements
ylabel('CO2 concentration (ppm)')
xlabel('time')
legend('Indoor', 'Indoor-controls','Outdoor')

% 
% numticks = get(gca,'XTick');
% dateticks = datenum(datenum(led.t.label)+numticks/86400);
% datestrings = datestr(dateticks,dateFormat);
% xticklabels(datestrings);
% 
% figure(3)
% plot(led.x.co2Air)
% hold on
% plot(led.d.co2Out)
% ylabel('CO2 concentration (mg m^{-3})')
% hold off
% legend('Indoor','Outdoor')
% 
% numticks = get(gca,'XTick');
% dateticks = datenum(datenum(led.t.label)+numticks/86400);
% datestrings = datestr(dateticks,dateFormat);
% xticklabels(datestrings);
% 

% 
% subplot(2,3,5)
% plot(led.d.iGlob)
% hold on
% plot(led.a.rParGhSun+led.a.rParGhLamp)
% plot(led.a.qLampIn)
% plot(led.a.rParGhSun)
% plot(led.a.rParGhLamp)
% legend('Outdoor ledobal solar radiation','PAR above the canopy (sun+lamp)',...
% 'Lamp electric input','PAR above the canopy (sun)', 'PAR above the canopy (lamp)')
% ylabel('W m^{-2}')
% 
% numticks = get(gca,'XTick');
% dateticks = datenum(datenum(led.t.label)+numticks/86400);
% datestrings = datestr(dateticks,dateFormat);
% xticklabels(datestrings);
% 
% subplot(2,3,6)
% plot(led.p.parJtoUmolSun*led.a.rParGhSun)
% hold on
% plot(led.p.zetaLampPar*led.a.rParGhLamp)
% legend('PPFD from the sun','PPFD from the lamp')
% ylabel('umol (PAR) m^{-2} s^{-1}')
% 
% numticks = get(gca,'XTick');
% dateticks = datenum(datenum(led.t.label)+numticks/86400);
% datestrings = datestr(dateticks,dateFormat);
% xticklabels(datestrings);

% Artificial Weather dataset
function weather = weatherDataset(length)
% make an artificial dataset to use as input for a createGreenLightModel
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

    length = ceil(length);
    weather = nan(length*288,9);
    time = 0:300:(length*86400-1);
    weather(:,1) = time;
    weather(:,2) = 350*max(0,sin(time*2*pi/86400));
    weather(:,3) = 5*sin(time*2*pi/86400)+15;
    weather(:,4) = 0.006*ones(length*288,1);
    weather(:,5) = co2ppm2dens(weather(:,3), 410);
    % weather(:,6) = 1*ones(length*288,1);
    weather(:,6) = 0;
    weather(:,7) = weather(:,3) - 20;
    weather(:,8) = 20*ones(length*288,1);

    % convert timestamps to datenum
    weather(:,1) = time/86400+1;
    weather(:,9) = dayLightSum(weather(:,1), weather(:,2));

end

%% Controls parameters
%  length - length of desired dataset(days)
%  controls will be matrix with 9 columns, in the following formats:
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
function controls = controlsParameters(length)
    length = ceil(length);
    controls = nan(length*288,9);
    time = 0:300:(length*86400-1);
    controls(:,1) = time;                           % timestamps
    controls(:,2) = 0;                              % 'thScr' Energy screen closure 			0-1 (1 is fully closed)
    controls(:,3) = 0;                              % 'blScr' Black out screen closure			0-1 (1 is fully closed)
    controls(:,4) = 1;                              % 'roof'  Average roof ventilation aperture	(average between lee side and wind side)	0-1 (1 is fully open)
    controls(:,5) = 5*sin(time*2*pi/86400)+25; %25;      % 'tPipe' Pipe rail temperature 			°C
    controls(:,6) = 5*sin(time*2*pi/86400)+25; %25;      % 'tGroPipe' Grow pipes temperature 			°C
    controls(:,7) = 1;                              % 'lamp' Toplights on/off                  0/1 (1 is on)
    controls(:,8) = 0;                              % 'intLamp' Interlight on/off                 0/1 (1 is on)
    controls(:,9) = 0;                              % 'extCo2' CO2 injection   

    % convert timestamps to datenum
    controls(:,1) = time/86400+1;
end

%% Indoor dataset
% make an artificial dataset to use as input for a createGreenLightModel
%   indoor          (optional) A 3 column matrix with:
%       indoor(:,1)     timestamps of the input [s] in regular intervals of 300, starting with 0
%       indoor(:,2)     temperature       [°C]             indoor air temperature
%       indoor(:,3)     vapor pressure    [Pa]             indoor vapor concentration
%       indoor(:,4)     co2 concentration [mg m^{-3}]      indoor co2 concentration
function indoor = indoorDataset(length)
    length = ceil(length);
    indoor = nan(length*288,4);
    time = 0:300:(length*86400-1);
    
    temperature_outside = 5*sin(time*2*pi/86400)-15;

    indoor(:,1) = time; 
    indoor(:,2) = 5*sin(time*2*pi/86400);  
    
    % Calculate saturation vapor pressure
    es = 6.112 * exp(17.67 * indoor(:,2) ./ (indoor(:,2) + 243.5));
    
    % Calculate indoor vapor pressure using RH
    indoor(:,3) = es * (70 / 100);
    
    indoor(:,4) = co2ppm2dens(temperature_outside, 410);

    % convert timestamps to datenum
    indoor(:,1) = time/86400+1;
end
