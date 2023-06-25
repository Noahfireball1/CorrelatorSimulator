classdef Navigation
    %Navigation defines properties associated with user-defined navigation parameters

    properties
        navFreq
        
    end

    properties (Access = private)
        sampleFreq
        simTime
    end

    properties (Dependent)
        navFreqSamples
        numUpdates
    end


    methods
        function obj = Navigation(config,simulation)

            obj.navFreq = config.frequency;
            obj.sampleFreq = simulation.sampleFreq;
            obj.simTime = simulation.time;
            

        end

        function navFreqSamples = get.navFreqSamples(obj)
            navFreqSamples = obj.sampleFreq/obj.navFreq;
        end
        function numUpdates = get.numUpdates(obj)
            numUpdates = obj.navFreq*obj.simTime;
        end
    end
end