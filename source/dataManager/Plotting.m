classdef Plotting
    %PLOTTING Summary of this class goes here
    %   Detailed explanation goes here

    properties
        docked = [];
        tiled = [];
    end

    methods

        function obj = Plotting(config)
            obj.docked = config.docked;
            obj.tiled = config.tiled;

        end

    end
end