classdef PLL2 < LoopFilters

    methods
        function obj = PLL2(pdiTime)
            % Assign PLL Parameters
            obj.assignProperties(bandwidth=2,coeffA=sqrt(2),order=2);

            % Calculate Gains 1,2
            obj.gain1 = obj.naturalFreq^2*(pdiTime) + obj.coeffA*obj.naturalFreq;

            obj.gain2 = -obj.coeffA*obj.naturalFreq;
        end

    end
end