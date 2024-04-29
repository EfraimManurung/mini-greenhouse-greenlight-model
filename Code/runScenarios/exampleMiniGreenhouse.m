% exampleMiniGreenhouse to run the GreenLight simulation
% Using createGreenLightModel

% Efraim Manurung, Information Technology Group
% Wageningen University
% efraim.efraimpartoginahotasi@wur.nl
% efraim.manurung@gmail.com

tic; % start the timer

% Weather argument for createGreenLightModel
seasonLength = 1; %5 it works when not using controls % season length in days
firstDay = 1; % days since beginning of data (01-01-2000)

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

[weather, startTime] = loadSelYearHiRes(firstDay, seasonLength);

secsInYear = seconds(startTime-datetime(year(startTime),1,1,0,0,0));

% number of seconds since beginning of year to startTime
weather(:,8) = soilTempNl(secsInYear+weather(:,1)); % add soil temperature

% Controls argument createGreenLightModel
using_controls = true;

if using_controls
    controls = makeArtificialControls(seasonLength);

    % Still using from existed database
    % time_controls = weather(:, 1);
    % time_controls = 0:300:86400;
    % air_temperature = weather(:, 3);
   
    % controls(:,1) = ;      % timestamps of the input [s] in regular intervals of 300, starting with 0
    % controls(:,2) = 0;                  % 'thScr' Energy screen closure 			0-1 (1 is fully closed)
    % controls(:,3) = 0;                  % 'blScr' Black out screen closure			0-1 (1 is fully closed)
    % controls(:,4) = 0;                  % 'roof'  Average roof ventilation aperture	(average between lee side and wind side)	0-1 (1 is fully open)
    % controls(:,5) = 5*sin(time*2*pi/86400)+15;    % 'tPipe' Pipe rail temperature 			°C
    % controls(:,6) = 5*sin(time*2*pi/86400)+15;    % 'tGroPipe' Grow pipes temperature 			°C
    % controls(:,7) = 1;                  % 'lamp' Toplights on/off                  0/1 (1 is on)
    % controls(:,8) = 0;                  % 'intLamp' Interlight on/off                 0/1 (1 is on)
    % controls(:,9) = 0;                  % 'extCo2' CO2 injection                     0/1 (1 is on)
end

%% Create an instance of createGreenLight with the default Vanthoor parameters
% led = createGreenLightModel('led', weather, startTime);
led = createGreenLightModel('led', weather, startTime, controls);

% Parameters for mini-greenhouse
% Parameters are good 
setParamsMiniGreenhouse(led);      % set greenhouse structure
setMiniGreenhouseLedParams(led);   % set lamp params

% Parameters for Bleiswjik2010, to make some tests
% setParamsBleiswijk2010(led);    % set greenhouse structure
% setBleiswijk2010LedParams(led); % set lamp params

% Set initial values for crop
led.x.cLeaf.val = 0.7*6240*10;
led.x.cStem.val = 0.25*6240*10;
led.x.cFruit.val = 0.05*6240*10;

%% Run simulation
solveFromFile(led, 'ode15s');

% set data to a fixed step size (5 minutes)
led = changeRes(led, 300);

%% Plot some outputs 
% see setGlAux, setGlStates, setGlInput to see more options

dateFormat = 'HH:00'; 
% This format can be changed, see help file for MATLAB function datestr

subplot(3,4,1)
plot(led.x.tAir)
hold on
plot(led.d.tOut)
ylabel('Temperature (°C)')
legend('Indoor','Outdoor')

numticks = get(gca,'XTick');
dateticks = datenum(datenum(led.t.label)+numticks/86400);
datestrings = datestr(dateticks,dateFormat);
xticklabels(datestrings);

subplot(3,4,2)
plot(led.x.vpAir)
hold on
plot(led.d.vpOut)
ylabel('Vapor pressure (Pa)')
legend('Indoor','Outdoor')

% numticks = get(gca,'XTick');
% dateticks = datenum(datenum(led.t.label)+numticks/86400);
% datestrings = datestr(dateticks,dateFormat);
% xticklabels(datestrings);

subplot(3,4,3)
plot(led.a.rhIn)
hold on
plot(100*vp2dens(led.d.tOut,led.d.vpOut)./rh2vaporDens(led.d.tOut,100));
ylabel('Relative humidity (%)')
legend('Indoor','Outdoor')

% numticks = get(gca,'XTick');
% dateticks = datenum(datenum(led.t.label)+numticks/86400);
% datestrings = datestr(dateticks,dateFormat);
% xticklabels(datestrings);

subplot(3,4,4)
plot(led.x.co2Air)
hold on
plot(led.d.co2Out)
ylabel('CO2 concentration (mg m^{-3})')
legend('Indoor','Outdoor')

% numticks = get(gca,'XTick');
% dateticks = datenum(datenum(led.t.label)+numticks/86400);
% datestrings = datestr(dateticks,dateFormat);
% xticklabels(datestrings);

subplot(3,4,5)
plot(led.a.co2InPpm)
hold on
plot(co2dens2ppm(led.d.tOut,1e-6*led.d.co2Out))
ylabel('CO2 concentration (ppm)')
legend('Indoor','Outdoor')

% numticks = get(gca,'XTick');
% dateticks = datenum(datenum(led.t.label)+numticks/86400);
% datestrings = datestr(dateticks,dateFormat);
% xticklabels(datestrings);

subplot(3,4,6)
plot(led.d.iGlob)
hold on
plot(led.a.rParGhSun+led.a.rParGhLamp)
plot(led.a.qLampIn)
plot(led.a.rParGhSun)
plot(led.a.rParGhLamp)
legend('Outdoor ledobal solar radiation','PAR above the canopy (sun+lamp)',...
'Lamp electric input','PAR above the canopy (sun)', 'PAR above the canopy (lamp)')
ylabel('W m^{-2}')

% numticks = get(gca,'XTick');
% dateticks = datenum(datenum(led.t.label)+numticks/86400);
% datestrings = datestr(dateticks,dateFormat);
% xticklabels(datestrings);

subplot(3,4,7)
plot(led.p.parJtoUmolSun*led.a.rParGhSun)
hold on
plot(led.p.zetaLampPar*led.a.rParGhLamp)
legend('PPFD from the sun','PPFD from the lamp')
ylabel('umol (PAR) m^{-2} s^{-1}')

% numticks = get(gca,'XTick');
% dateticks = datenum(datenum(led.t.label)+numticks/86400);
% datestrings = datestr(dateticks,dateFormat);
% xticklabels(datestrings);

subplot(3,4,9)
plot(led.a.mcAirCan)
hold on
plot(led.a.mcAirBuf)
plot(led.a.mcBufAir)
plot(led.a.mcOrgAir)
ylabel('mg m^{-2} s^{-1}')
legend('Net assimilation (CO_2)', 'Net photosynthesis (gross photosynthesis minus photorespirattion, CH_2O)',...
'Growth respiration (CH_2O)','Maintenance respiration (CH_2O)')

% numticks = get(gca,'XTick');
% dateticks = datenum(datenum(led.t.label)+numticks/86400);
% datestrings = datestr(dateticks,dateFormat);
% xticklabels(datestrings);

subplot(3,4,10)
plot(led.x.cFruit)
hold on
plot(led.x.cStem)
plot(led.x.cLeaf)
plot(led.x.cBuf)
ylabel('mg (CH_2O) m^{-2}')
yyaxis right
plot(led.a.lai)
ylabel('m^2 m^{-2}')
legend('Fruit dry weight','Stem dry weight','Leaf dry weight','Buffer content','LAI')

% numticks = get(gca,'XTick');
% dateticks = datenum(datenum(led.t.label)+numticks/86400);
% datestrings = datestr(dateticks,dateFormat);
% xticklabels(datestrings);

subplot(3,4,11)
plot(led.x.cFruit)
ylabel('mg (CH_2O) m^{-2}')
yyaxis right
plot(led.a.mcFruitHar)
ylabel('mg (CH_2O) m^{-2} s^{-1}')
legend('Fruit dry weight','Fruit harvest')

% numticks = get(gca,'XTick');
% dateticks = datenum(datenum(led.t.label)+numticks/86400);
% datestrings = datestr(dateticks,dateFormat);
% xticklabels(datestrings);

toc;

%% make an artificiual dataset to use as input for a createGreenLightModel
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
function controls = makeArtificialControls(length)
    length = ceil(length);
    controls = nan(length*288,9);
    time = 0:300:(length*86400-1);
    controls(:,1) = time;                           % timestamps
    controls(:,2) = 0;                              % 'thScr' Energy screen closure 			0-1 (1 is fully closed)
    controls(:,3) = 0;                              % 'blScr' Black out screen closure			0-1 (1 is fully closed)
    controls(:,4) = 0;                              % 'roof'  Average roof ventilation aperture	(average between lee side and wind side)	0-1 (1 is fully open)
    controls(:,5) = 5*sin(time*2*pi/86400)+15;      % 'tPipe' Pipe rail temperature 			°C
    controls(:,6) = 5*sin(time*2*pi/86400)+15;      % 'tGroPipe' Grow pipes temperature 			°C
    controls(:,7) = 1;                              % 'lamp' Toplights on/off                  0/1 (1 is on)
    controls(:,8) = 0;                              % 'intLamp' Interlight on/off                 0/1 (1 is on)
    controls(:,9) = 0;                              % 'extCo2' CO2 injection    
end


