classdef NavigationDataClass < handle

    properties
        position
        velocity
        clockBias
        clockDrift
        rangeResolution
        dopplerResolution
        receiveTime
        cno
        losPos
        losVel
        losClockBias
        losClockDrift
        covariance
        rangeResolutionVariance
        dopplerResolutionVariance
    end

    methods
        function obj = NavigationDataClass(channelNum,dataLength)
            obj.position = NaN(3,dataLength);
            obj.velocity = NaN(3,dataLength);
            obj.clockBias = NaN(1,dataLength);
            obj.clockDrift = NaN(1,dataLength);
            obj.rangeResolution = NaN(channelNum,dataLength);
            obj.dopplerResolution = NaN(channelNum,dataLength);
            obj.receiveTime = NaN(1,dataLength);
            obj.cno = NaN(channelNum,dataLength);
            obj.losPos = NaN(3,dataLength);
            obj.losVel = NaN(3,dataLength);
            obj.losClockBias = NaN(1,dataLength);
            obj.losClockDrift = NaN(1,dataLength);
            obj.covariance = NaN(8,dataLength);
            obj.rangeResolutionVariance = NaN(channelNum,dataLength);
            obj.dopplerResolutionVariance = NaN(channelNum,dataLength);
        end

    end
end