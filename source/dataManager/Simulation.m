classdef Simulation
    %SIMULATION defines settings specific to the correlator simulation

    properties

        sampleFreq
        numSeconds

    end

    properties (Constant)
        C = 299792458;
        omege_ie = 7.2921151467e-5;
        intermedFreq = 4.092e6;
        chipFreq = 1.023e6;
        gpsL1Freq = 1575.42e6
        gpsPRNArray = 1:32;
        pdiTimeCodePeriods = 20;
        pdiTime = 20e-3;
        codeOffset = 0.5;
        codeLength = 1023;

    end

    properties (Dependent)

        dataLength
        numSamples
        gpsL1WaveLength
        gpsPRNSequences
        chipWidth

    end

    methods
        function obj = Simulation(config)

            obj.sampleFreq = config.sampleFreq;
            obj.numSeconds = config.time;

        end

        function dataLength = get.dataLength(obj)
            dataLength = obj.numSamples/obj.pdiTime;
        end
        function numSamples = get.numSamples(obj)
            numSamples = obj.numSeconds*obj.sampleFreq;
        end
        function gpsL1WaveLength = get.gpsL1WaveLength(obj)
            gpsL1WaveLength = obj.C/obj.gpsL1Freq;
        end
        function chipWidth = get.chipWidth(obj)
            chipWidth = obj.C/obj.chipFreq;
        end
        function gpsPRNSequences = get.gpsPRNSequences(obj)
            G1 = ones(1,10);
            G2 = ones(1,10);

            G1out = zeros(1,1023);
            G2out = zeros(1,1023);
            gpsPRNSequences = nan([37,obj.codeLength]);

            for prn = 1:37
                for i = 1:obj.codeLength
                    G1sum = mod(G1(10)+G1(3),2);
                    switch prn
                        case 1
                            S1 = G2(2);
                            S2 = G2(6);
                        case 2
                            S1 = G2(3);
                            S2 = G2(7);
                        case 3
                            S1 = G2(4);
                            S2 = G2(8);
                        case 4
                            S1 = G2(5);
                            S2 = G2(9);
                        case 5
                            S1 = G2(1);
                            S2 = G2(9);
                        case 6
                            S1 = G2(2);
                            S2 = G2(10);
                        case 7
                            S1 = G2(1);
                            S2 = G2(8);
                        case 8
                            S1 = G2(2);
                            S2 = G2(9);
                        case 9
                            S1 = G2(3);
                            S2 = G2(10);
                        case 10
                            S1 = G2(2);
                            S2 = G2(3);
                        case 11
                            S1 = G2(3);
                            S2 = G2(4);
                        case 12
                            S1 = G2(5);
                            S2 = G2(6);
                        case 13
                            S1 = G2(6);
                            S2 = G2(7);
                        case 14
                            S1 = G2(7);
                            S2 = G2(8);
                        case 15
                            S1 = G2(8);
                            S2 = G2(9);
                        case 16
                            S1 = G2(9);
                            S2 = G2(10);
                        case 17
                            S1 = G2(1);
                            S2 = G2(4);
                        case 18
                            S1 = G2(2);
                            S2 = G2(5);
                        case 19
                            S1 = G2(3);
                            S2 = G2(6);
                        case 20
                            S1 = G2(4);
                            S2 = G2(7);
                        case 21
                            S1 = G2(5);
                            S2 = G2(8);
                        case 22
                            S1 = G2(6);
                            S2 = G2(9);
                        case 23
                            S1 = G2(1);
                            S2 = G2(3);
                        case 24
                            S1 = G2(4);
                            S2 = G2(6);
                        case 25
                            S1 = G2(5);
                            S2 = G2(7);
                        case 26
                            S1 = G2(6);
                            S2 = G2(8);
                        case 27
                            S1 = G2(7);
                            S2 = G2(9);
                        case 28
                            S1 = G2(8);
                            S2 = G2(10);
                        case 29
                            S1 = G2(1);
                            S2 = G2(6);
                        case 30
                            S1 = G2(2);
                            S2 = G2(7);
                        case 31
                            S1 = G2(3);
                            S2 = G2(8);
                        case 32
                            S1 = G2(4);
                            S2 = G2(9);
                        case 33
                            S1 = G2(5);
                            S2 = G2(10);
                        case 34
                            S1 = G2(4);
                            S2 = G2(10);
                        case 35
                            S1 = G2(1);
                            S2 = G2(7);
                        case 36
                            S1 = G2(2);
                            S2 = G2(8);
                        case 37
                            S1 = G2(4);
                            S2 = G2(10);
                        otherwise
                            fprintf('Specified PRN Number not found. Please try again. \n')
                    end
                    phaseSelectorOut = mod(S1+S2,2);
                    G2sum = mod(G2(10)+G2(9)+G2(8)+G2(6)+G2(3)+G2(2),2);

                    G1out(i) = G1(10);
                    G2out(i) = G2(10);

                    gpsPRNSequences(prn,i) = mod(phaseSelectorOut + G1out(i),2);

                    G1 = [G1sum G1(1:9)];
                    G2 = [G2sum G2(1:9)];

                end
                zero2one = gpsPRNSequences(prn,:) == 0;
                gpsPRNSequences(prn,zero2one) = -1;
            end
        end
    end
end