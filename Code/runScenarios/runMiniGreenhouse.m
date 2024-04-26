% Author: Efraim Manurung
% Information Technology Group, Wageningen University
%
% efraim.efraimpartoginahotasi@wur.nl 
% efraim.manurung@gmail.com
%
% Based on the David Katzin paper
% David Katzin, Simon van Mourik, Frank Kempkes, and Eldert J. Van Henten. 2020. 
% “GreenLight - An Open Source Model for Greenhouses with Supplemental Lighting: 
% Evaluation of Heat Requirements under LED and HPS Lamps.” Biosystems Engineering 194: 61–81. 
% https://doi.org/10.1016/j.biosystemseng.2020.03.010
%
% Modifying GreenLight
% Using runGreenLight is simple, but it doesn't allow to modify the GreenLight model. 
% In order to modify the model before simulating it, a
% more complex approach is needed. This can be done in the following way:

% 7.3.31. Setting up a GreenLight instance using createGreenLightModel.m
% The function createGreenLightModel creates a DynamicModel object that
% represents a greenhouse. THe function requires the following arguments:

% - lampType :  same as in runGreenLight.m (Section 7.3.1)
% - weather  :  same as in runGreenLight.m (Section 7.3.1), with the 
%               exception that the timestamps (column 1) are given in seconds since the 
%               beginning of the data (starting in 0), and that elevation (column 10) is not provided.
% - startTime: a time label (in MATLAB%s datenum format) indicating the date and time when the data started.

function gl = runMiniGreenhouse(lampType, weather, filename, paramNames, paramVals, isMature, controls, indoor)
% runMiniGreenhouse Create an instance of the GreenLight model and run a simulation. 
% Uses parameter settings that represent a SETBLEIWSWIJK2010LEDPARAMS
% Modify parameters for a GreenLight.
% (see also setBleiswijk2010LedParams(gl)).
% Inputs:
%   gl - a DynamicModel object to be use as a GreenLight model
%
% Inputs (all optional):
%   lampType    Lamp type to simulate, may be 'hps', 'led', or 'none'. 
%               Default value is 'none'     
%   weather     Weather inputs for the model. If this argument empty, artificial 
%               weather data for a 5 day seasons will be generated. If this
%               argument is a scalar number, artifical data will be
%               generated for as many days as this scalar number.
%               Otherwise, weather needs to be a 9-column matrix 
%               in the following format:
%                   weather(:,1)    timestamps of the input [datenum] in 5 minute intervals
%                   weather(:,2)    radiation     [W m^{-2}]  outdoor global irradiation 
%                   weather(:,3)    temperature   [°C]        outdoor air temperature
%                   weather(:,4)    humidity      [kg m^{-3}] outdoor vapor concentration
%                   weather(:,5)    co2 [kg{CO2} m^{-3}{air}] outdoor CO2 concentration
%                   weather(:,6)    wind        [m s^{-1}] outdoor wind speed
%                   weather(:,7)    sky temperature [°C]
%                   weather(:,8)    temperature of external soil layer [°C]
%                   weather(:,9)    daily radiation sum [MJ m^{-2} day^{-1}]
%                   weather(:,10)   elevation [m above sea level] (optional, default is 0)
%   filename    Name of file where the simulation results will be saved.
%               If empty or blank (''), no file will be saved.
%   paramNames  Array of strings with names of parameters to modify beyond
%               their default values. Example: ["lampsOn" "lampsOff"]
%               NOTE: "dependent parameters" (those defined in
%               setDepParams) should be changed by changing their defining
%               parameters. For example setting aPipe as 0 should be done
%               by setting lPipe or phiPipeE as 0. See setDepParams for a
%               list of dependent parameters.
%   paramVals   Array of values (corresponding to paramNames) with modified
%               parameter values. Example: [2 16]
%   isMature    If true, simulation will start with a mature crop. Default
%               is false.
%   controls    (optional) A matrix with 8 columns, in the following format:
%                   controls(:,1)     timestamps of the input [s] in regular intervals of 300, starting with 0
%                   controls(:,2)     Energy screen closure 			0-1 (1 is fully closed)
%                   controls(:,3)     Black out screen closure			0-1 (1 is fully closed)
%                   controls(:,4)     Average roof ventilation aperture	(average between lee side and wind side)	0-1 (1 is fully open)
%                   controls(:,5)     Pipe rail temperature 			°C
%                   controls(:,6)     Grow pipes temperature 			°C
%                   controls(:,7)     Toplights on/off                  0/1 (1 is on)
%                   controls(:,8)     Interlight on/off                 0/1 (1 is on)
%                   controls(:,9)     CO2 injection                     0/1 (1 is on)
%
%   indoor      (optional) A 3 column matrix with:
%                   indoor(:,1)     timestamps of the input [s] in regular intervals of 300, starting with 0
%                   indoor(:,2)     temperature       [°C]             indoor air temperature
%                   indoor(:,3)     vapor pressure    [Pa]             indoor vapor concentration
%                   indoor(:,4)     co2 concentration [mg m^{-3}]      indoor co2 concentration   
%
%   Output:
%   gl          An instance of the GreenLight model with the completed
%               simulation. Data of this simulation is given in 5 minute
%               intervals.

    %% Set default values
    absTol = 1e-6; % default is 1e-6
    relTol = 1e-3; % default 

    % For lampType
    if exist('lampType', 'var')
        if strcmpi(lampType, 'hps')
            lampType = 'hps';
        elseif strcmpi(lampType, 'led')
            lampType = 'led';
        else
            lampType = 'none';
        end
    else
        lampType = 'none';
    end
    
    % For weather
    if ~exist('weather', 'var') || isempty(weather)
        weather = 5;    
    end
    
    if isscalar(weather)
        weather = makeArtificialWeather(weather);
    end

    if size(weather,2) < 10
        elevation = 0;
    else
        elevation = weather(1,10);
    end

    % For filename
    if ~exist('filename', 'var') || isempty(filename)
        filename = '';
    end
    if ~exist('mature', 'var') || isempty(isMature)
        isMature = false;
    end
    if strcmp(filename,'')
        saveFile = false;
    else
        saveFile = true;
    end

    % For controls
    if ~exist('controls', 'var') || isempty(controls)
        controls = 5;    
    end
    
    if isscalar(controls)
        controls = makeArtificialControls(controls);
        % controls = ............(controls);
    end
    
    % For indoor
    if ~exist('indoor', 'var') || isempty(indoor)
        indoor = 5;    
    end
    
    if isscalar(indoor)
        % indoor = ............(indoor);
        indoor = makeArtificialIndoor(indoor);
    end
    
    disp(datetime('now'))
    disp(['starting runGreenLight. Lamp type: ' lampType '. Filename: ' filename '.']);
    tic;
    
    %% Prepare and run simulation
    % convert weather timestamps from datenum to seconds since beginning of
    % data
    startTime = datetime(weather(1,1),'ConvertFrom', 'datenum');
    weather(:,1) = (weather(:,1)-weather(1,1))*86400;

    gl = createGreenLightModel(lampType, weather, startTime, controls, indoor);
    setParamsBleiswijk2010(gl);     % set 
    setBleiswijk2010LedParams(gl);  % set lamp params
    setParam(gl, 'hElevation', elevation); % set elevation

    % Modify parameters based on function input
    % ParamName should be a vector of strings
    % paramVal should be a respective numeric vector
    if exist('paramNames', 'var') && ~isempty(paramNames)
        for k=1:length(paramNames)
            setParam(gl, paramNames(k), paramVals(k));
        end
    end

    % Reset other parameters that depend on previously defined parameters
    setDepParams(gl);

    %% Start with a mature crop
    if isMature
        setVal(gl.x.cFruit, 2.8e5);
        setVal(gl.x.cLeaf, 0.9e5);
        setVal(gl.x.cStem, 2.5e5);
        setVal(gl.x.tCanSum, 3000);
    end

    % Run simulation
    options = odeset('AbsTol', absTol, 'RelTol', relTol);
    solveFromFile(gl, 'ode15s', options);
    gl = changeRes(gl, 300);

    if saveFile
        save(filename);
    end

    toc;

end

function indoor = makeArtificialIndoor(length)
%   indoor      (optional) A 3 column matrix with:
%                   indoor(:,1)     timestamps of the input [s] in regular intervals of 300, starting with 0
%                   indoor(:,2)     temperature       [°C]             indoor air temperature
%                   indoor(:,3)     vapor pressure    [Pa]             indoor vapor concentration
%                   indoor(:,4)     co2 concentration [mg m^{-3}]      indoor co2 concentration   
    time = 0:300:(length*86400-1);
    indoor(:,1) = time;
    indoor(:,2) = 20*ones(length*288,1);
    indoor(:,3) = 1000;
    indoor(:,4) = 400;
end

function controls = makeArtificialControls(length)
% make an artifical dataset to use as input for a GreenLight instance
%   length  - length of desired dataset (days)
%         controls(:,1)     timestamps of the input [s] in regular intervals of 300, starting with 0
%         controls(:,2)     Energy screen closure 			0-1 (1 is fully closed)
%         controls(:,3)     Black out screen closure			0-1 (1 is fully closed)
%         controls(:,4)     Average roof ventilation aperture	(average between lee side and wind side)	0-1 (1 is fully open)
%         controls(:,5)     Pipe rail temperature 			°C
%         controls(:,6)     Grow pipes temperature 			°C
%         controls(:,7)     Toplights on/off                  0/1 (1 is on)
%         controls(:,8)     Interlight on/off                 0/1 (1 is on)
%         controls(:,9)     CO2 injection                     0/1 (1 is on)
    time = 0:300:(length*86400-1);
    controls(:,1) = time;
    controls(:,2) = 0;
    controls(:,3) = 0;
    controls(:,4) = 0;
    controls(:,5) = 25.0;
    controls(:,6) = 25.0;
    controls(:,7) = 0;
    controls(:,8) = 0;
    controls(:,9) = 0;

     % convert timestamps to datenum
    controls(:,1) = time/86400+1;
end

function weather = makeArtificialWeather(length)
% make an artifical dataset to use as input for a GreenLight instance
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
    weather(:,6) = 1*ones(length*288,1);
    weather(:,7) = weather(:,3) - 20;
    weather(:,8) = 20*ones(length*288,1);
    
    % convert timestamps to datenum
    weather(:,1) = time/86400+1;
    weather(:,9) = dayLightSum(weather(:,1), weather(:,2));

end