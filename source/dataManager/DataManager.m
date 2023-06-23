classdef DataManager < handle
    %DATAMANAGER Parses, Loads, and Creates settings needed to run the correlator sim

    properties (Access = public)
        directories         = UserDirectories();
        general             = General();
        simulation          = Simulation();
        trajectory          = Trajectory();
        processing          = Processing();
        plotting            = Plotting();
        constants           = Constants();

    end

    methods (Access = public)
        function obj = DataManager(config)

        end

    end

    methods(Access = private)

    end
end