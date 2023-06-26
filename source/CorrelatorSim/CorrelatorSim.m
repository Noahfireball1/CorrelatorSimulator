classdef CorrelatorSim < handle
    %CORRELATORSIM simulates on the correlator level and then performs scalar or vector tracking

    properties (Access = public)

    end
    properties (Access = public)
        dir
        gen
        sim
        traj
        nav
        plotting
    end

    methods
        function obj = CorrelatorSim(settings)
            obj.grabSettings(settings)

        end

        function calcSVPos(obj)

            eph = obj.gen.ephemeris;
            transmitTime = obj.sim.initialTransmitTime;
            transitTime = 0.068;

            SatellitePositions(eph,transmitTime,transitTime)
            
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