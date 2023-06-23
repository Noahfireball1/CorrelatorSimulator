classdef Simulation
    %SIMULATION Summary of this class goes here
    %   Detailed explanation goes here

    properties

        dataLength
        numSamples
        
    end

    properties (Constant)

        sampleFreq
        intermedFreq
        chipFreq
        gpsL1Freq
        gpsL1Lambda
        gpsPRNArray
        pdiTimeCodePeriods
        pdiTime
        codeOffset
        chipWidth
        codeLength

    end

    methods
        function obj = Simulation(inputArg1,inputArg2)
            %SIMULATION Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end