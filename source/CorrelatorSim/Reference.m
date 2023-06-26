classdef Reference < dynamicprops
    %REFERENCE Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
        
    end

    methods
        function obj = Reference(sim)

            numSats = size(sim.satellitePositions.svPosX,1);

            for sv = 1:numSats

                % Add a channel
                channelName = sprintf('channel%i',sv);
                obj.(channelName) = obj.addprop(channelName);

                % Initialize Channel
                obj.(channelName) = Channel(sim,sv);
                
                
            end

        end

    end
end