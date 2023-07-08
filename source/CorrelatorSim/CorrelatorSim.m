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
    properties (Dependent)
        numChannels
        dataLength
    end

    methods (Access = public)
        function obj = CorrelatorSim(settings)
            obj.grabSettings(settings)

            % Calculating Initial Satellite Positions
            obj.calcSVPos;
            % Initialize Reference Data
            obj.initializeReference;
            % initialize Estimated Data
            obj.initializeEstimate;
            % Initialize Loop Filters
            obj.initializeLoopFilters;
            % Initialize Scalar Data
            obj.initializeScalar;
            % Initialize Navigation Data
            obj.initializeNavigation;

        end

        function step(obj)
            fileEnd = 0;
            count = 1;
            firstLoop = 1;
            sampleCount = 0;
            navCount = 1;

            textprogressbar('[correlator-sim] Simulation Progress: ')

            while fileEnd == 0

                % Processing next sample
                if mod(count,10) == 0
                    try
                        textprogressbar(count/10)
                    catch
                        textprogressbar('[correlator-sim] Simulation Progress: ')
                    end
                end
                
                % Kalman Time Update
                [~,obj.reference] = TimeUpdate(obj.reference,1/obj.sim.sampleFreq,count);

                % Calculate Correlators
                PropagateCorrelators(obj,sampleCount)

                % Perform Navigation Solution


                if sampleCount > obj.sim.dataLength
                    fileEnd = 1;
                end
                count = count + 1;
                if count > 1000
                    fileEnd = 1;
                end
            end
            fprintf('\n')
            textprogressbar('[correlator-sim] Simulation Complete!')
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
        function calcSVPos(obj)

            eph = obj.gen.ephemeris;
            transmitTime = obj.sim.initialTransmitTime;
            position = obj.traj.position;

            obj.satellitePositions = SatellitePositions(eph,transmitTime,position);

        end

        function initializeReference(obj)

            navLength = obj.nav.numUpdates;

            obj.reference = Reference(obj,obj.numChannels,obj.dataLength,navLength);
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

            obj.scalar = Scalar(obj.numChannels,obj.dataLength);

        end

        function initializeNavigation(obj)

            obj.navigation = NavigationDataClass(obj.numChannels,obj.dataLength);

        end

    end

    methods
        function numChannels = get.numChannels(obj)
            numChannels = size(obj.satellitePositions.svPosX,1);
        end
        function dataLength = get.dataLength(obj)
            dataLength = obj.sim.dataLength;
        end
    end
end