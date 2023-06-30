classdef LoopFilters < handle
    %LOOPFILTERS Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = protected)
        bandWidth
        naturalFreq
        coeffA
        coeffB
    end

    properties (Access = public)
        gain1
        gain2
    end

    methods (Access = protected)

        function assignProperties(obj,nameValue)
            arguments
                obj
                nameValue.coeffA
                nameValue.coeffB
                nameValue.bandWidth
                nameValue.order
            end

            obj.bandWidth = nameValue.bandWidth;
            obj.coeffA = nameValue.coeffA;
            if nameValue.order == 3
                obj.coeffB = nameValue.coeffB;
                obj.naturalFreq = obj.bandWidth/0.7845;
            else
                obj.naturalFreq = obj.bandWidth/0.53;
            end

        end

    end
end