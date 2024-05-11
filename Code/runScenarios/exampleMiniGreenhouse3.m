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
v.rhAir = DynamicElement('v.rhAir', [floor(indoor(:,1)) indoor(:,3)]);
v.co2Air = DynamicElement('v.co2Air', [floor(indoor(:,1)) indoor(:,4)]);

% number of seconds since beginning of year to startTime
secsInYear = seconds(startTime-datetime(year(startTime),1,1,0,0,0));

%outdoor(:,7) = skyTempRdam(outdoor(:,3), datenum(startTime)+outdoor(:,1)/86400); % add sky temperature
outdoor(:,7) = outdoor(:,3) - 10;
outdoor(:,8) = soilTempNl(secsInYear+outdoor(:,1)); % add soil temperature

%% Create an instance of createGreenLight with the default Vanthoor parameters
% led = createGreenLightModel('led', outdoor, startTime, controls, indoor);
led = createGreenLightModel(lampType, outdoor, startTime);
led_controls = createGreenLightModel(lampType, outdoor, startTime,controls);

% Parameters for mini-greenhouse
setParamsMiniGreenhouse(led);      % set greenhouse structure
setMiniGreenhouseLedParams(led);   % set lamp params

setParamsMiniGreenhouse(led_controls);      % set greenhouse structure
setMiniGreenhouseLedParams(led_controls);   % set lamp params
%% Control parameters
% Read setGIParams.m about the explanation and default values of the control parameters
% setParam(led, 'rhMax', 87);        % upper bound on relative humidity  
% setParam(led_controls, 'rhMax', 87);        % upper bound on relative humidity  

% Set initial values for crop
% start with 3.12 plants/m2, assume they are each 2 g = 6240 mg/m2.
% Check the setGlinit.m for more information
% Default values
led.x.cLeaf.val = 0.7*6240;
led.x.cStem.val = 0.25*6240;
led.x.cFruit.val = 0.05*6240;

% Set initial values for crop
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
rrmseRhair = sqrt(mean((led_controls.a.rhIn.val(1:compareLength,2)-v.rhAir.val(1:compareLength,2)).^2))./mean(v.rhAir.val(1:compareLength,2));
rrmseCo2air  = sqrt(mean((led_controls.x.co2Air.val(1:compareLength,2)-v.co2Air.val(1:compareLength,2)).^2))./mean(v.co2Air.val(:,2));

rmseTair = sqrt(mean((led_controls.x.tAir.val(1:compareLength,2) - v.tAir.val(:,2)).^2));
rmseRhair = sqrt(mean((led_controls.a.rhIn.val(1:compareLength,2)-v.rhAir.val(1:compareLength,2)).^2));
rmseCo2air = sqrt(mean((led_controls.x.co2Air.val(1:compareLength,2) - v.co2Air.val(1:compareLength,2)).^2));

% Calculate ME 
meTair = mean(led_controls.x.tAir.val(1:compareLength,2) - v.tAir.val(:,2));
meRhair = mean(led_controls.a.rhIn.val(1:compareLength,2)-v.rhAir.val(1:compareLength,2));
meCo2air = mean(led_controls.x.co2Air.val(1:compareLength,2) - v.co2Air.val(1:compareLength,2));

disp('rrmseTair [%]'); disp(rrmseTair);
disp('rrmseRhair[%]'); disp(rrmseRhair);
disp('rrmseCo2air [%]'); disp(rrmseCo2air);

disp('rmseTair [°C]'); disp(rmseTair);
disp('rmseRhair [%]'); disp(rmseRhair);
disp('rmseCo2air [ppm]'); disp(rmseCo2air);

disp('meTair [°C]'); disp(meTair);
disp('meRhair[%]'); disp(meRhair);
disp('meCo2air [ppm]'); disp(meCo2air);

%% Plot some outputs 
% see setGlAux, setGlStates, setGlInput to see more options

% dateFormat = 'HH:00'; 
% This format can be changed, see help file for MATLAB function datestr

%% Figure 1 TEMPERATURE
figure(1)
plot(led.x.tAir,'LineWidth',1.5)
hold on
plot(led_controls.x.tAir,'LineWidth',1.5)
plot(v.tAir.val(:,1), v.tAir.val(:,2),'LineWidth',1.5)
plot(led_controls.d.tOut,'LineWidth',1.5)
hold off
xlabel('Time')
ylabel('Temperature (°C)')
legend('Indoor','Indoor-controls','Indoor-real measurements','Outdoor')
setXAxisTicksAndLabels(led.t.label, seasonLength)

%% Figure 2 RELATIVE HUMIDITY
figure(2)
plot(led.a.rhIn,'LineWidth',1.5)
hold on
plot(led_controls.a.rhIn, 'LineWidth',1.5)
plot(v.rhAir.val(:,1), v.rhAir.val(:,2),'LineWidth',1.5)
plot(100*vp2dens(led.d.tOut,led.d.vpOut)./rh2vaporDens(led.d.tOut,100),'LineWidth',1.5);
hold off
xlabel('Time')
ylabel('Relative humidity (%)')
legend('Indoor', 'Indoor-controls','Indoor-real measurements','Outdoor')
setXAxisTicksAndLabels(led.t.label, seasonLength)

%% Figure 3 CO2 IN PPM
figure(3)
plot(led.a.co2InPpm,'LineWidth',1.5)
hold on
plot(led_controls.a.co2InPpm,'LineWidth',1.5)
plot(v.co2Air.val(:,1), v.co2Air.val(:,2),'LineWidth',1.5)
plot(co2dens2ppm(led.d.tOut,1e-6*led.d.co2Out),'LineWidth',1.5)
hold off
xlabel('Time')
ylabel('CO2 concentration (ppm)')
legend('Indoor', 'Indoor-controls','Indoor-real measurements','Outdoor')
setXAxisTicksAndLabels(led.t.label, seasonLength)

%% Figure 4 PPFD
figure(4)
plot(led.p.parJtoUmolSun * led.a.rParGhSun,'LineWidth',1.5)
hold on
plot(led.p.zetaLampPar * led.a.rParGhLamp,'LineWidth',1.5)
plot(led_controls.p.zetaLampPar * led_controls.a.rParGhLamp,'LineWidth',1.5)
hold off
xlabel('Time')
ylabel('umol (PAR) m^{-2} s^{-1}')
legend('PPFD from the sun', 'PPFD from the lamp', 'PPFD from the lamp - controls')
setXAxisTicksAndLabels(led.t.label, seasonLength)

%% FUnction for the figures
function setXAxisTicksAndLabels(timeLabels, seasonLength)
    numTicks = get(gca,'XTick');
    dateticks = datenum(datenum(timeLabels) + numTicks / (60*60*24)); % Assuming timestamps are in seconds
    if seasonLength > 2
        datestrings = datestr(dateticks, 'dd');
    else
        datestrings = datestr(dateticks, 'HH:00');
    end
    xticks(numTicks);
    xticklabels(datestrings);
    xtickangle(45);
end