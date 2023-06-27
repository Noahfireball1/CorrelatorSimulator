classdef EstimateChannel < handle
    %CHANNEL Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
        prn
        cnoL1
        losVel
        carrFreq
        doppler
        codeFreq
        codePhase
        phase
        codeError
        oldCarrierError
        oldCodeError
        carrierError
        oldFreqError
        oldCodeFreq
        oldCarrierFreq
        receiveTime
        transmitTime
        updateCounter
        codePeriodCounter
        initialReceiveTime
        initialTransmitTime
        dopplerRate
        oldCarrier
        newCodePeriodStarted
    end
    methods (Access = public)
        function obj = EstimateChannel(sim,sv)

        end
    end

    methods (Access = private)

        function transmitTime = calcTransmitTime(obj,svProps,userTraj,sv)
            range = vecnorm([svProps.svPosX(sv) svProps.svPosY(sv) svProps.svPosZ(sv)] - userTraj.position,2,2);
            signalEffects = ((range + userTraj.ionosphereDelay + userTraj.troposphereDelay + userTraj.clockBias - svProps.svClockCorr(sv)*svProps.C)/svProps.C);

            transmitTime = obj.receiveTime - signalEffects;
        end

        function velocity = calcVelocity(~,userTraj,svProps,sv)
            svVelVec = [svProps.svVelX(sv) svProps.svVelY(sv) svProps.svVelZ(sv)];
            userVel = userTraj.velocity;

            velocity = svProps.unitVectors(sv,:)*(svVelVec - userVel)';
        end

        function carrierFrequency = calcCarrierFreq(obj,svProps,userTraj,freqs)
            dynamicsFreq = (1 - obj.losVel/svProps.C);
            driftFreq = (1 + userTraj.clockDrift);

            carrierFrequency = dynamicsFreq*freqs.gpsL1Freq/driftFreq - freqs.gpsL1Freq + freqs.intermedFreq;
        end

        function codePhase = calcCodePhase(obj,freqs)
            partialCodePeriodSec = mod(obj.transmitTime,1e-3);

            codePhase = partialCodePeriodSec*freqs.chipFreq;
        end

        function phase = calcPhase(~,svProps,userTraj,freqs,sv)
            range = vecnorm([svProps.svPosX(sv) svProps.svPosY(sv) svProps.svPosZ(sv)] - userTraj.position,2,2);

            phase = -(range + userTraj.ionosphereDelay + userTraj.troposphereDelay + userTraj.clockBias - svProps.svClockCorr(sv)*svProps.C)/freqs.gpsL1WaveLength;
        end

    end
end