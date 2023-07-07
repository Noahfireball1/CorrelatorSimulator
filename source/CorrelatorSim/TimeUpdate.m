classdef TimeUpdate < handle

    properties (Access = public)
        phi
        processNoise = zeros(8,8)
    end

    properties (Access = private)
        C = 299792458;
        sigmaP = 2;
        h0 = 2.8e-23;
        h2 = 2.4e-23;
    end

    methods (Access = public)
        function [obj,reference] = TimeUpdate(reference,timeStep,count)

            obj.calcProcess(timeStep);
            
            reference.stateVector = obj.phi*reference.stateVector;
            reference.position(:,count) = reference.stateVector(1:3);
            reference.velocity(:,count) = reference.stateVector(4:6);
            reference.clockBias(count) = reference.stateVector(7);
            reference.clockDrift(count) = reference.stateVector(8);

        end

        function calcProcess(obj,timeStep)
            obj.phi = [1 0 0 timeStep 0 0 0 0; ...
                0 1 0 0 timeStep 0 0 0; ...
                0 0 1 0 0 timeStep 0 0; ...
                0 0 0 1 0 0 0 0;        ...
                0 0 0 0 1 0 0 0;        ...
                0 0 0 0 0 1 0 0;        ...
                0 0 0 0 0 0 1 timeStep; ...
                0 0 0 0 0 0 0 1];

            a = timeStep^3/3;
            b = timeStep^2/2;
            c = timeStep^2/2;
            d = timeStep;

            varB = obj.C^2*obj.h0/2;
            varD = obj.C^2*2*pi^2*obj.h2;

            obj.processNoise(1,1) = obj.sigmaP*a;
            obj.processNoise(1,4) = obj.sigmaP*b;
            obj.processNoise(2,2) = obj.sigmaP*a;
            obj.processNoise(2,5) = obj.sigmaP*b;
            obj.processNoise(3,3) = obj.sigmaP*a;
            obj.processNoise(3,6) = obj.sigmaP*b;
            obj.processNoise(4,4) = obj.sigmaP*d;
            obj.processNoise(4,1) = obj.sigmaP*c;
            obj.processNoise(5,5) = obj.sigmaP*d;
            obj.processNoise(5,2) = obj.sigmaP*c;
            obj.processNoise(6,6) = obj.sigmaP*d;
            obj.processNoise(6,3) = obj.sigmaP*c;
            obj.processNoise(7:8,7:8) = varD*[a b;c d];
            obj.processNoise(7,7) = obj.processNoise(7,7) + varB*timeStep;

        end
    end

end