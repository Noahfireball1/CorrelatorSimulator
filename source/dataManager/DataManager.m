classdef DataManager < handle
    %DATAMANAGER Parses, Loads, and Creates settings needed to run the correlator sim

    properties (Access = public)
        directories         = [];
        general             = [];
        simulation          = [];
        trajectory          = [];
        navigation          = [];
        plotting            = [];

    end

    methods (Access = public)
        function obj = DataManager(configFilePath)
            [~,configName,~] = fileparts(configFilePath);

            config = ReadYaml(configFilePath);

            obj.directories = UserDirectories();
            if ~exist(append(obj.directories.data,configName),"dir")
                mkdir(append(obj.directories.data,configName))
            end
            obj.directories.data = append(obj.directories.data,configName,filesep);

            obj.general = General(config.general,obj.directories);
            obj.simulation = Simulation(config.simulation,obj.general);
            obj.trajectory = Trajectory(config.trajectory);
            obj.navigation = Navigation(config.navigation,config.simulation);
            obj.plotting = Plotting(config.plotting);

        end

    end

    methods(Access = private)

    end
end