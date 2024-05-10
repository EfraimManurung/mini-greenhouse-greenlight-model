% exampleMiniGreenhouse3 to run the GreenLight simulation
% Using createGreenLightModel
%
% Efraim Manurung, Information Technology Group
% Wageningen University
% efraim.efraimpartoginahotasi@wur.nl
% efraim.manurung@gmail.com
%
% Based on:
% David Katzin, Simon van Mourik, Frank Kempkes, and Eldert J. Van Henten. 2020. 
% "GreenLight - An Open Source Model for Greenhouses with Supplemental Lighting: Evaluation of Heat Requirements under LED and HPS Lamps.” 
% Biosystems Engineering 194: 61–81. https://doi.org/10.1016/j.biosystemseng.2020.03.010

tic; % start the timer
%% Set up the model
% Weather argument for createGreenLightModel
seasonLength = 1; % season length in days
firstDay = 1; % days since beginning of data 

% Choice of lamp
lampType = 'led'; % 'led', 'hps', or 'none'

[outdoor, indoor, controls, startTime] = loadMiniGreenhouseData2(firstDay, seasonLength);

% DynamicElements for the measured data
v.tAir = DynamicElement('v.tAir', [floor(indoor(:,1)) indoor(:,2)]);
%v.vpAir = DynamicElement('v.vpAir', [floor(indoor(:,1)) indoor(:,3)]);
v.co2Air = DynamicElement('v.co2Air', [floor(indoor(:,1)) indoor(:,4)]);

secsInYear = seconds(startTime-datetime(year(startTime),1,1,0,0,0));
% number of seconds since beginning of year to startTime

outdoor(:,7) = skyTempRdam(outdoor(:,3), datenum(startTime)+outdoor(:,1)/86400); % add sky temperature
outdoor(:,8) = soilTempNl(secsInYear+weather(:,1)); % add soil temperature

%% Create an instance of createGreenLight with the default Vanthoor parameters
% convert weather timestamps from datenum to seconds since beginning of data
% startTime = datetime(weather(1,1),'ConvertFrom','datenum');

% led = createGreenLightModel('led', weather, startTime, controls, indoor);
led = createGreenLightModel(lampType, weather, startTime);
led_controls = createGreenLightModel(lampType, weather, startTime,controls);

% Parameters for mini-greenhouse
setParamsMiniGreenhouse(led);      % set greenhouse structure
setMiniGreenhouseLedParams(led);   % set lamp params

setParamsMiniGreenhouse(led_controls);      % set greenhouse structure
setMiniGreenhouseLedParams(led_controls);   % set lamp params
%% Control parameters
% Read setGIParams.m about the explanation and default values of the control parameters
% setParam(led, 'rhMax', 60);        % upper bound on relative humidity  

% Set initial values for crop
% start with 3.12 plants/m2, assume they are each 2 g = 6240 mg/m2.
% Check the setGlinit.m for more information
% Default values
led.x.cLeaf.val = 0.7*6240;
led.x.cStem.val = 0.25*6240;
led.x.cFruit.val = 0.05*6240;

% % Set initial values for crop
led_controls.x.cLeaf.val = 0.7*6240;
led_controls.x.cStem.val = 0.25*6240;
led_controls.x.cFruit.val = 0.05*6240;

%% Run simulation
solveFromFile(led, 'ode15s');
solveFromFile(led_controls, 'ode15s');

% set data to a fixed step size (5 minutes)
led = changeRes(led, 300);
led_controls = changeRes(led_controls, 300);

toc;

%% Get RRMSEs between simulation and measurements
% Check that the measured data and the simulations have the same size. 
% If one of them is bigger, some data points of the longer dataset will be
% discarded
mesLength = length(v.tAir.val(:,1)); % the length (array size) of the measurement data
simLength = length(led_controls.x.tAir.val(:,1)); % the length (array size) of the simulated data
compareLength = min(mesLength, simLength);

rrmseTair = sqrt(mean((led_controls.x.tAir.val(1:compareLength,2)-v.tAir.val(:,2)).^2))./mean(v.tAir.val(1:compareLength,2));
%rrmseVpair = sqrt(mean((led.x.vpAir.val(1:compareLength,2)-v.vpAir.val(1:compareLength,2)).^2))./mean(v.vpAir.val(1:compareLength,2));
rrmseCo2air  = sqrt(mean((led_controls.x.co2Air.val(1:compareLength,2)-v.co2Air.val(1:compareLength,2)).^2))./mean(v.co2Air.val(:,2));


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
%plot(led_controls.d.tOut,'LineWidth',1.5)
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
%plot(100*vp2dens(led_controls.d.tOut,led_controls.d.vpOut)./rh2vaporDens(led_controls.d.tOut,100),'LineWidth',1.5);
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
%plot(co2dens2ppm(led_controls.d.tOut,1e-6*led_controls.d.co2Out),'LineWidth',1.5)
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


%% Figure 4 PPFD
figure(4)
plot(led.p.parJtoUmolSun*led.a.rParGhSun,'LineWidth',1.5)
hold on
% plot(led_controls.p.parJtoUmolSun*led_controls.a.rParGhSun)
plot(led.p.zetaLampPar*led.a.rParGhLamp,'LineWidth',1.5)
plot(led_controls.p.zetaLampPar*led_controls.a.rParGhLamp,'LineWidth',1.5)
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
ylabel('umol (PAR) m^{-2} s^{-1}')
xlabel('time')
legend('PPFD from the sun','PPFD from the lamp', 'PPFD from the lamp - controls')