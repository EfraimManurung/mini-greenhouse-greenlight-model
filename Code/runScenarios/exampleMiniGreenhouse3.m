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
seasonLength = 5; % season length in days
firstDay = 1; % days since beginning of data 

% Choice of lamp
lampType = 'led'; 

[outdoor, indoor, controls, startTime] = loadMiniGreenhouseData2(firstDay, seasonLength);

% DynamicElements for the measured data
v.tAir = DynamicElement('v.tAir', [floor(indoor(:,1)) indoor(:,2)]);
v.rhAir = DynamicElement('v.rhAir', [floor(indoor(:,1)) indoor(:,3)]);
v.co2Air = DynamicElement('v.co2Air', [floor(indoor(:,1)) indoor(:,4)]);
v.iInside = DynamicElement('v.iInside', [floor(indoor(:,1)) indoor(:,5)]);

% number of seconds since beginning of year to startTime
secsInYear = seconds(startTime-datetime(year(startTime),1,1,0,0,0));

%outdoor(:,7) = skyTempRdam(outdoor(:,3), datenum(startTime)+outdoor(:,1)/86400); % add sky temperature
outdoor(:,7) = outdoor(:,3) - 10;
outdoor(:,8) = soilTempNl(secsInYear+outdoor(:,1)); % add soil temperature

%% Create an instance of createGreenLight with the default Vanthoor parameters
% led = createGreenLightModel('led', outdoor, startTime, controls, indoor);
led = createGreenLightModel(lampType, outdoor, startTime, controls);

% Parameters for mini-greenhouse
setParamsMiniGreenhouse(led);      % set greenhouse structure
setMiniGreenhouseLedParams(led);   % set lamp params
%% Control parameters
% Read setGIParams.m about the explanation and default values of the control parameters
% setParam(led, 'rhMax', 50);        % upper bound on relative humidity  

% Set initial values for crop
% start with 3.12 plants/m2, assume they are each 2 g = 6240 mg/m2.
% Check the setGlinit.m for more information
% Default values    
led.x.cLeaf.val = 0.7*6240;     
led.x.cStem.val = 0.25*6240;    
led.x.cFruit.val = 0.05*6240;   

% Default values
% led.x.cLeaf.val = 0.01;     
% led.x.cStem.val = 0.01;    
% led.x.cFruit.val = 0.01; 

%% Run simulation
solveFromFile(led, 'ode15s');

% set data to a fixed step size (5 minutes)
led = changeRes(led, 300);

toc;
%% Get RRMSEs between simulation and measurements
% Check that the measured data and the simulations have the same size. 
% If one of them is bigger, some data points of the longer dataset will be
% discarded
mesLength = length(v.tAir.val(:,1)); % the length (array size) of the measurement data
simLength = length(led.x.tAir.val(:,1)); % the length (array size) of the simulated data
compareLength = min(mesLength, simLength);

% Apply the multiplier to led.a.rhIn values
multiplier = 0.61; %0.85; %0.61; %0.83;
if exist('multiplier', 'var') && ~isempty(multiplier)
    led.a.rhIn.val(:,2) = led.a.rhIn.val(:,2) * multiplier;
end

% Added PAR from sun and lamp
sunLampIrradiance = (led.a.rParGhSun.val(1:compareLength,2)+led.a.rParGhLamp.val(1:compareLength,2)) - 1.66;

% Calculate RRMSE
rrmseTair = (sqrt(mean((led.x.tAir.val(1:compareLength,2)-v.tAir.val(:,2)).^2))./mean(v.tAir.val(1:compareLength,2))) * 100;
rrmseRhair = (sqrt(mean((led.a.rhIn.val(1:compareLength,2)-v.rhAir.val(1:compareLength,2)).^2))./mean(v.rhAir.val(1:compareLength,2))) * 100;
rrmseCo2air  = (sqrt(mean((led.x.co2Air.val(1:compareLength,2)-v.co2Air.val(1:compareLength,2)).^2))./mean(v.co2Air.val(1:compareLength,2))) * 100;
rrmseIinside = (sqrt(mean((sunLampIrradiance - v.iInside.val(1:compareLength,2)).^2))./mean(v.iInside.val(1:compareLength,2))) * 100;

% Calculate RMSE
rmseTair = sqrt(mean((led.x.tAir.val(1:compareLength,2) - v.tAir.val(1:compareLength,2)).^2));
rmseRhair = sqrt(mean((led.a.rhIn.val(1:compareLength,2)-v.rhAir.val(1:compareLength,2)).^2));
rmseCo2air = sqrt(mean((led.x.co2Air.val(1:compareLength,2) - v.co2Air.val(1:compareLength,2)).^2));
rmseIinside = sqrt(mean((sunLampIrradiance - v.iInside.val(1:compareLength,2)).^2));

% Calculate ME 
meTair = mean(led.x.tAir.val(1:compareLength,2) - v.tAir.val(:,2));
meRhair = mean(led.a.rhIn.val(1:compareLength,2)- v.rhAir.val(1:compareLength,2));
meCo2air = mean(led.x.co2Air.val(1:compareLength,2) - v.co2Air.val(1:compareLength,2));
meIinside = mean(sunLampIrradiance - v.iInside.val(1:compareLength,2));

% Save the output 
save exampleMiniGreenhouse

% Display the values
fprintf('\n');
if exist('multiplier', 'var') && ~isempty(multiplier)
    fprintf('Multiplier RH: %.2f\n', multiplier);
end
fprintf('Season Length: %d day(s) \n', seasonLength);
fprintf('---------------------------------------------\n');
fprintf('| Metric          | Value       | Unit       \n');
fprintf('---------------------------------------------\n');
fprintf('| RRMSE Tair      | %-12.2f| %%              \n', rrmseTair);
fprintf('| RRMSE Rhair     | %-12.2f| %%              \n', rrmseRhair);
fprintf('| RRMSE Co2air    | %-12.2f| %%              \n', rrmseCo2air);
fprintf('| RRMSE IInside   | %-12.2f| %%              \n', rrmseIinside);
fprintf('| RMSE Tair       | %-12.2f| °C              \n', rmseTair);
fprintf('| RMSE Rhair      | %-12.2f| %%              \n', rmseRhair);
fprintf('| RMSE Co2air     | %-12.2f| ppm             \n', rmseCo2air);
fprintf('| RMSE IInside    | %-12.2f| W m^{-2}        \n', rmseIinside);
fprintf('| ME Tair         | %-12.2f| °C              \n', meTair);
fprintf('| ME Rhair        | %-12.2f| %%              \n', meRhair);
fprintf('| ME Co2air       | %-12.2f| ppm             \n', meCo2air);
fprintf('| ME Iinside      | %-12.2f| W m^{-2}        \n', meIinside);
fprintf('---------------------------------------------\n');

%% Plot some outputs 
% see setGlAux, setGlStates, setGlInput to see more options
% Create a new figure
figure;

%% TEMPERATURE FIGURES
subplot(4, 1, 1); % 3 rows, 1 column, 1st subplot
plot(v.tAir.val(:, 1), v.tAir.val(:, 2), 'LineWidth', 1.0);
hold on;
plot(led.x.tAir, 'LineWidth', 1.0);
hold off;
ylabel('Air Temp In [°C]');
legend('Indoor-Measured', 'Indoor-Simulated');
setXAxisTicksAndLabels(led.t.label, seasonLength);

subplot(4, 1, 2);
plot(led.d.tOut, 'LineWidth', 1.0);
ylabel('Air Temp Out [°C]');
legend('Outdoor-Measured');
setXAxisTicksAndLabels(led.t.label, seasonLength);

subplot(4, 1, 3);
plot(led.d.iGlob, 'LineWidth', 1.0); 
ylabel('Illumin Out [W m^{-2}]');
legend('Outdoor-Measured');
setXAxisTicksAndLabels(led.t.label, seasonLength);

subplot(4, 1, 4);
plot(controls(:,1), controls(:,4), 'LineWidth', 1.0);
ylabel('Fans [-]');
ylim([-0.05 1.05]);
legend('Controls-Measured');
setXAxisTicksAndLabels(led.t.label, seasonLength);

%% RELATIVE HUMIDITY FIGURES
figure;
subplot(4, 1, 1); % 3 rows, 1 column, 2nd subplot
plot(v.rhAir.val(:, 1), v.rhAir.val(:, 2), 'LineWidth', 1.0);
hold on;
plot(led.a.rhIn, 'LineWidth', 1.0);
hold off;
ylabel('Rel Humid In [%]');
legend('Indoor-Measured', 'Indoor-Simulated');
setXAxisTicksAndLabels(led.t.label, seasonLength);

subplot(4, 1, 2);
plot(100 * vp2dens(led.d.tOut, led.d.vpOut) ./ rh2vaporDens(led.d.tOut, 100), 'LineWidth', 1.0); 
ylabel('Rel Humid Out [%]');
legend('Outdoor-Measured');
setXAxisTicksAndLabels(led.t.label, seasonLength);

subplot(4, 1, 3); 
plot(v.tAir.val(:, 1), v.tAir.val(:, 2), 'LineWidth', 1.0);
ylabel('Air Temp In [°C]');
legend('Indoor-Measured');
setXAxisTicksAndLabels(led.t.label, seasonLength);

subplot(4, 1, 4);
plot(controls(:,1), controls(:,4), 'LineWidth', 1.0);
ylabel('Fans [-]');
ylim([-0.05 1.05]);
legend('Controls-Measured');
setXAxisTicksAndLabels(led.t.label, seasonLength);

%% PAR light level FIGURES
figure;

subplot(3, 1, 1);
plot(v.iInside.val(:, 1), v.iInside.val(:, 2), 'LineWidth', 1.0);
hold on;
plot(led.a.rParGhSun + led.a.rParGhLamp, 'LineWidth', 1.0);
hold off;
ylabel('PAR in [W m^{-2}]');
legend('Measured - PAR above the canopy (sun+lamp)', 'Simulated - PAR above the canopy (sun+lamp)');
setXAxisTicksAndLabels(led.t.label, seasonLength);

% Create a subplot and plot with green leaf color
subplot(3, 1, 2);
% plot(led.d.iGlob, 'LineWidth', 1.0, 'Color', [0.13, 0.55, 0.13]); % RGB values for green leaf color
plot(led.d.iGlob, 'LineWidth', 1.0) 
ylabel('Illumin Out [W m^{-2}]');
legend('Outdoor-Measured');
setXAxisTicksAndLabels(led.t.label, seasonLength);

subplot(3, 1, 3);
plot(controls(:,1), controls(:,7), 'LineWidth', 1.0);
ylabel('Lamps [-]');
ylim([-0.05 1.05]);
legend('Controls-Measured');
setXAxisTicksAndLabels(led.t.label, seasonLength);


%% Clear the workspace
clear;

%% Function for the figures
function setXAxisTicksAndLabels(timeLabels, seasonLength)
    numTicks = get(gca,'XTick');
    dateticks = datenum(datenum(timeLabels) + numTicks / (60*60*24)); % Assuming timestamps are in seconds
    if seasonLength > 2
        datestrings = datestr(dateticks, 'dd');
        xlabel('Time [d]')
    else
        xlabel('Time [h]')
        datestrings = datestr(dateticks, 'HH:00');
    end
    xticks(numTicks);
    xticklabels(datestrings);
    %xtickangle(45);
end

