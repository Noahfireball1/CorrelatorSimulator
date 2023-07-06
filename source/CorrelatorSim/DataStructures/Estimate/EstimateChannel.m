classdef EstimateChannel < handle
    %CHANNEL Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
        prn
        cno
        losVel
        carrierFreq
        doppler
        codeFreq
        codePhase
        phase
        codeError = 0;
        oldCarrierError = 0;
        oldCodeError = 0;
        carrierError = 0;
        oldFreqError = 0;
        oldCodeFreq
        oldCarrierFreq
        receiveTime
        transmitTime
        updateCounter = 1;
        codePeriodCounter = 0;
        initialReceiveTime
        initialCodePhase
        initialTransmitTime
        dopplerRate
        oldCarrier
        newCodePeriodStarted = 0;
        filter
    end

    properties (Access = private)
        powerL1 = 45;
    end


    methods (Access = public)
        function obj = EstimateChannel(sim,phaseError,sv)
            svProps = sim.satellitePositions;
            userTraj = sim.traj;
            simFreq = sim.sim;
            ref = sim.reference;

            obj.prn = svProps.ID(sv);
            obj.cno = 10^(obj.powerL1/10);
            obj.losVel = obj.calcVelocity(userTraj,svProps,sv);
            obj.carrierFreq = obj.calcCarrierFreq(svProps,userTraj,simFreq);
            obj.doppler = obj.carrierFreq - simFreq.intermedFreq;
            obj.codeFreq = obj.doppler*simFreq.chipFreq/simFreq.gpsL1Freq + simFreq.chipFreq;
            obj.codePhase = ref.channel1.codePhase - phaseError;
            obj.phase = ref.channel1.phase - obj.calcInitialCarrPhaseError(phaseError,sim);
            obj.oldCodeFreq = obj.codeFreq;
            obj.oldCarrierFreq = obj.carrierFreq;
            obj.receiveTime = ref.channel1.receiveTime;
            obj.transmitTime =  ref.channel1.transmitTime - mod(ref.channel1.transmitTime,1e-3);
            obj.initialReceiveTime = obj.receiveTime;
            obj.initialCodePhase = obj.codePhase;
            obj.initialTransmitTime = obj.transmitTime;
            obj.dopplerRate = -obj.phase;
            obj.oldCarrier = obj.dopplerRate;
            obj.filter = EstimateFilter(sim,obj.carrierFreq);



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

        function initialCarrPhaseError = calcInitialCarrPhaseError(~,phaseError,sim)
            initialRangeError = phaseError/sim.sim.chipFreq*sim.sim.C;
            initialCarrPhaseError = mod((initialRangeError*sim.reference.channel1.carrierFreq/sim.sim.C*2*pi),2*pi);
            if initialCarrPhaseError > pi
                initialCarrPhaseError = initialCarrPhaseError - 2*pi;
            end
        end


    end

    methods

    end
end