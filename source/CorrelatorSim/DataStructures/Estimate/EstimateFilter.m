classdef EstimateFilter < handle

    properties (Access = public)
        dynamics % A
        inputs   % B
        observability % C
        state % X
        kalmanGain = [0.541080410673989;7.83430131462034;56.6176707167131] % L
    end

    methods (Access = public)
        function obj = EstimateFilter(sim,carrierFreq)
            obj.dynamics = [1, sim.sim.pdiTime, sim.sim.pdiTime^2/2; 0, 1, sim.sim.pdiTime; 0, 0, 1];
            obj.inputs = [-sim.sim.pdiTime;0;0];
            obj.observability = [1, sim.sim.pdiTime/2, sim.sim.pdiTime^2/6];
            obj.state = [0;carrierFreq;0];
        end
    end
end