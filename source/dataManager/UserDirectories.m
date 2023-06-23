classdef UserDirectories
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        project = [];
        config = [];
        data = [];
        output = [];
        source = [];
    end

    methods
        function obj = UserDirectories()

            obj.project = fileparts(which("run_simulation.m"));
            obj.config = append(obj.project,filesep,'config',filesep);
            obj.data = append(obj.project,filesep,'data',filesep);
            obj.output = append(obj.project,filesep,'output',filesep);
            obj.source = append(obj.project,filesep,'source',filesep);

        end
    end
end