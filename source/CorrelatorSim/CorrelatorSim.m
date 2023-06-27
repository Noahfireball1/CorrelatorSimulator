classdef CorrelatorSim < handle
    %CORRELATORSIM simulates on the correlator level and then performs scalar or vector tracking

    properties (Access = public)
        satellitePositions
        reference
        estimate

    end
    properties (Access = public)
        dir
        gen
        sim
        traj
        nav
        plotting
    end

    methods (Access = public)
        function obj = CorrelatorSim(settings)
            obj.grabSettings(settings)

        end

        function calcSVPos(obj)

            eph = obj.gen.ephemeris;
            transmitTime = obj.sim.initialTransmitTime;
            position = obj.traj.position;

            obj.satellitePositions = SatellitePositions(eph,transmitTime,position);
            
        end

        function initializeReference(obj)

            obj.reference = Reference(obj);
        end

        function initializeEstimate(obj)

            obj.estimate = Estimate(obj);
        end
    end

    methods (Access = private)

        function grabSettings(obj,settings)

            obj.dir = settings.directories;
            obj.gen = settings.general;
            obj.sim = settings.simulation;
            obj.traj = settings.trajectory;
            obj.nav = settings.navigation;
            obj.plotting = settings.plotting;

        end

    end
end