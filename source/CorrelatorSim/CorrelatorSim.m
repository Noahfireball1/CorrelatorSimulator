classdef CorrelatorSim < handle
    %CORRELATORSIM simulates on the correlator level and then performs scalar or vector tracking

    properties (Access = public)
        satellitePositions
        reference
        estimate
        loopFilters
        scalar
        navigation

    end
    properties (Access = public)
        dir
        gen
        sim
        traj
        nav
        plotting
    end
    properties (Access = private)
        numChannels
        dataLength
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

        function initializeLoopFilters(obj)
            obj.loopFilters = LoopFiltersDataClass();

            obj.loopFilters.PLL3 = PLL3(obj.sim.pdiTime);
            obj.loopFilters.PLL2 = PLL2(obj.sim.pdiTime);
            obj.loopFilters.DLL3 = DLL3(obj.sim.pdiTime);
            obj.loopFilters.DLL2 = DLL2(obj.sim.pdiTime);
            obj.loopFilters.FLL2 = FLL2(obj.sim.pdiTime);
        end

        function initializeScalar(obj)

            obj.scalar = Scalar();

        end

        function initializeNavigation(obj)

            obj.navigation = NavigationDataClass();

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