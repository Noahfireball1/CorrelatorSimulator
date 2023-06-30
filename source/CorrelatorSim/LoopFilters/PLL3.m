classdef PLL3 < LoopFilters

    properties (Access = public)
        gain3
    end

    methods
        function obj = PLL3(pdiTime)

            % Assign PLL Parameters
            obj.assignProperties(bandwidth=10,coeffA=1.1,coeffB=2.4,order=3);

            % Calculate Gains 1,2,3
            obj.gain1 = obj.naturalFreq^3*(pdiTime^2) + obj.coeffA*obj.naturalFreq^2*pdiTime + obj.coeffB*obj.naturalFreq;

            obj.gain2 = -(obj.coeffA*obj.naturalFreq^2*(pdiTime) + 2*obj.coeffB*obj.naturalFreq);

            obj.gain3 = obj.coeffB*obj.naturalFreq;

        end
    end

end