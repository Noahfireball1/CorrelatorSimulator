classdef EstimateFilter < handle

    properties (Access = public)
        dynamics % A
        inputs   % B
        observability % C
        state % X
        kalmanGain % L
    end

    methods (Access = public)
        function obj = EstimateFilter(sim,sv)
        end
    end
end