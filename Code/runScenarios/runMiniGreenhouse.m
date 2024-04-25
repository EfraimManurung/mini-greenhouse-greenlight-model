% Author: Efraim Manurung
% Based on the David Katzin PhD Thesis
% 
% efraim.efraimpartoginahotasi@wur.nl 
% efraim.manurung@gmail.com
%
% Modifying GreenLight
% Using runGreenLight is simple, but it doesn't allow to modify the GreenLight model. 
% In order to modify the model before simulating it, a
% more complex approach is needed. This can be done in the following way:

% 7.3.31. Setting up a GreenLight instance using createGreenLightModel.m
% The function createGreenLightModel creates a DynamicModel object that
% represents a greenhouse. THe function requires the following arguments:

% - lampType: same as in runGreenLight.m (Section 7.3.1)

function gl = runMiniGreenhouse(lampType, weather, filename, control, paramNames, paramVals, isMature)
% runMiniGreenhouse Create an instance of the GreenLight model and run a simulation. 
% Uses parameter settings that represent a modern set lamp parameters according to modern specifications.
% (see also 
