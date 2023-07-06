classdef TimeUpdate < handle

    properties
        Property1
    end

    methods (Access = public)
        function [obj,reference] = TimeUpdate(reference,timeStep)
            obj.calcProcess(timeStep);


        end

        function calcProcess(obj,timeStep)


        end
    end
end