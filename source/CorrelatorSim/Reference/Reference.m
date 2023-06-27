classdef Reference < dynamicprops

    methods
        function obj = Reference(sim)

            numSats = size(sim.satellitePositions.svPosX,1);

            for sv = 1:numSats

                % Add a channel
                channelName = sprintf('channel%i',sv);
                obj.(channelName) = obj.addprop(channelName);

                % Initialize Channel
                obj.(channelName) = ReferenceChannel(sim,sv);
                
                
            end

        end

    end
end