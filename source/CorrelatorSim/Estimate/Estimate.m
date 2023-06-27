classdef Estimate < dynamicprops

    methods
        function obj = Estimate(sim)

            numSats = size(sim.satellitePositions.svPosX,1);

            for sv = 1:numSats

                % Add a Channel,Filter
                channelName = sprintf('channel%i',sv);
                filterName = sprintf('filter%i',sv);
                obj.(channelName) = obj.addprop(channelName);
                obj.(filterName) = obj.addprop(filterName);

                % Initialize Channel,Filter
                obj.(channelName) = EstimateChannel(sim,sv);
                obj.(filterName) = EstimateFilter(sim,sv);
                
                
            end

        end

    end
end