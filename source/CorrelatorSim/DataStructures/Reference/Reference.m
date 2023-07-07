classdef Reference < dynamicprops

    properties (Access = public)
        carrierPhase
        codePhase
        carrierFrequency
        codeFrequency
        receiveTime
        transmitTime

        elevation
        azimuth
        rangeError
        positionError
        position
        velocity
        clockBias
        clockDrift

        stateVector
    end

    methods
        function obj = Reference(sim,numChannels,dataLength,navLength)

            obj.carrierPhase = NaN(numChannels,dataLength);
            obj.codePhase = NaN(numChannels,dataLength);
            obj.carrierFrequency = NaN(numChannels,dataLength);
            obj.codeFrequency = NaN(numChannels,dataLength);
            obj.receiveTime = NaN(numChannels,dataLength);
            obj.transmitTime = NaN(numChannels,dataLength);

            obj.elevation = NaN(numChannels,navLength);
            obj.azimuth = NaN(numChannels,navLength);
            obj.rangeError = NaN(numChannels,navLength);
            obj.positionError = NaN(3,navLength);
            obj.position = NaN(3,navLength);
            obj.velocity = NaN(3,navLength);
            obj.clockBias = NaN(1,navLength);
            obj.clockDrift = NaN(1,navLength);

            obj.stateVector = [sim.traj.position';sim.traj.velocity';sim.traj.clockBias;sim.traj.clockDrift];

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