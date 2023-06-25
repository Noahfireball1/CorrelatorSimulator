%% Formatting
clc
clear
close all
format shortg

useGUI = 0;

%% Adding Paths Based on User's Directories
projectRoot = fileparts(which(mfilename));
addpath(genpath(projectRoot))
configDir = append(projectRoot,filesep,'config',filesep);
dataDir = append(projectRoot,filesep,'data',filesep);
outputDir = append(projectRoot,filesep,'output',filesep);

%% Select a Configuration File
if useGUI
    inputFile = uigetfile('*.yaml','Select Input File',configDir);
    inputFilePath = append(configDir,inputFile);
else
    inputFilePath = append(configDir,'example.yaml');
end

%% Initializing Simulation
settings = DataManager(inputFilePath);
% sim = CorrelatorSim(settings);