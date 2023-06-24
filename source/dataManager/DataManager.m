classdef DataManager < handle
    %DATAMANAGER Parses, Loads, and Creates settings needed to run the correlator sim

    properties (Access = public)
        directories         = [];
        general             = [];
        simulation          = [];
        trajectory          = [];
        navigation          = [];
        plotting            = [];
        constants           = Constants();

    end

    methods (Access = public)
        function obj = DataManager(configFilePath)

            config = ReadYaml(configFilePath);

            obj.directories = UserDirectories();
            obj.general = General(config.general,obj.directories);
            obj.simulation = Simulation(config.simulation);
            obj.trajectory = Trajectory(config.trajectory);
            obj.navigation = Navigation(config.navigation);
            obj.plotting = Plotting(config.plotting);

        end

    end

    methods(Access = private)

    end
end