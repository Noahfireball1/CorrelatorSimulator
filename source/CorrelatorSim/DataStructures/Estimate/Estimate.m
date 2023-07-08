classdef Estimate < dynamicprops

    properties(Access = public)
        position_ecef
        position_lla
        velocity_ecef
        clockBias = 0;
        clockDrift = 0;
        stateVector
        stateCovariance = diag(10000*ones(1,8));
        

    end

    properties (Access = private)
        initPosError = [0;0;0];
        initVeloError = [0;0;0];
        initRangeError;
    end

    methods
        function obj = Estimate(sim)

            obj.initRangeError = norm(obj.initPosError);

            obj.position_ecef = sim.traj.position' - obj.initPosError;
            obj.position_lla = ecef2lla(obj.position_ecef','WGS84');
            obj.velocity_ecef = sim.traj.velocity' + obj.initVeloError;
            obj.stateVector = [obj.position_ecef;obj.velocity_ecef;obj.clockBias;obj.clockDrift];

            numSats = size(sim.satellitePositions.svPosX,1);
            for sv = 1:numSats

                % Add a Channel,Filter
                channelName = sprintf('channel%i',sv);
                obj.(channelName) = obj.addprop(channelName);

                % Initialize Channel,Filter
                phaseError = obj.calcPhaseError(sim,sv);
                obj.(channelName) = EstimateChannel(sim,phaseError,sv);
                obj.(channelName).filter = EstimateFilter(sim,sv);

            end

        end

        function phaseError = calcPhaseError(obj,sim,sv)
            posErr = sim.traj.position' - obj.position_ecef;
            rangeErr = sim.satellitePositions.unitVectors(sv,:)*(posErr);
            totTimeErr = rangeErr/sim.sim.C;
            phaseError = totTimeErr*sim.sim.chipFreq;
        end

    end

end