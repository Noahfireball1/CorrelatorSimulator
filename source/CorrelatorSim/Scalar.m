classdef Scalar < handle

    properties
        carrierError
        averageCarrierError
        codeError
        carrierFrequencyError
        codeFrequencyError
        carrierDiscriminator
        codeDiscriminator
        IP
        QP
        IE
        QE
        IL
        QL
        receiveTime
        transmitTime
        rangeResolution

    end

    methods
        function obj = Scalar(channelNum,dataLength)
            obj.carrierError = NaN(channelNum,dataLength);
            obj.averageCarrierError = NaN(channelNum,dataLength);
            obj.codeError = NaN(channelNum,dataLength);
            obj.carrierFrequencyError = NaN(channelNum,dataLength);
            obj.codeFrequencyError = NaN(channelNum,dataLength);
            obj.carrierDiscriminator = NaN(channelNum,dataLength);
            obj.codeDiscriminator = NaN(channelNum,dataLength);
            obj.IP = NaN(channelNum,dataLength);
            obj.QP = NaN(channelNum,dataLength);
            obj.IE = NaN(channelNum,dataLength);
            obj.QE = NaN(channelNum,dataLength);
            obj.IL = NaN(channelNum,dataLength);
            obj.QL = NaN(channelNum,dataLength);
            obj.receiveTime = NaN(channelNum,dataLength);
            obj.transmitTime = NaN(channelNum,dataLength);
            obj.rangeResolution = NaN(channelNum,dataLength);
        end

    end
end