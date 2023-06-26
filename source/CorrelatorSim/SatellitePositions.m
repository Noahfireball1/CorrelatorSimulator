classdef SatellitePositions < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
       svPosX
       svPosY
       svPosZ
       svVelX
       svVelY
       svVelZ
       svClockCorr
       unitVectors
       transmitTime
       transitTime
       ID

    end

    properties (Access = private)
        TGD
        transmissionTime
        svClockDriftRate
        svClockDrift
        svClockBias
        Crc
        Crs
        Cuc
        Cus
        Cic
        Cis
        deltaN
        M0
        eccentricity
        sqrtA
        toe
        omega0
        i0
        omegaDot
        iDot
    end

    properties (Constant)
        gpsPi = 3.1415926535898;        % Pi used in the GPS coordinate system
        Omegae_dot = 7.2921151467e-5;   % Earth rotation rate, [rad/s]
        GM = 3.986005e14;               % Earth's universal [m^3/s^2]
        F = -4.442807633e-10;           % Constant, [sec/(meter)^(1/2)]
        C = 299792458                   % Speed of light [m/s]
    end

    methods
        function obj = SatellitePositions(eph,transmitTime,position)

            % Assign Properties from ephemeris file
            obj.associateEph(eph)

            % Initialize transitTime to be 68 ms for all satellites
            transitTime = 0.068*ones([length(obj.toe) 1]);

            for i = 1:10

                % Calculate Satellite Positions
                obj.calcSVPositions(transmitTime,transitTime);

                % Recalculate Transit Time
                svPos = [obj.svPosX obj.svPosY obj.svPosZ];
                transitTime = vecnorm((svPos - position),2,2)/obj.C;

            end

            % Calculate user-sv unit vectors
            obj.calcUnitVectors(position)

            obj.transmitTime = transmitTime;
            obj.transitTime = transitTime;


        end
    end
    methods (Access = private)

        function associateEph(obj,eph)
            obj.TGD = eph.TGD;
            obj.transmissionTime = eph.TransmissionTime;
            obj.svClockDriftRate = eph.SVClockDriftRate;
            obj.svClockDrift = eph.SVClockDrift;
            obj.svClockBias = eph.SVClockBias;
            obj.Crc = eph.Crc;
            obj.Crs = eph.Crs;
            obj.Cuc = eph.Cuc;
            obj.Cus = eph.Cus;
            obj.Cic = eph.Cic;
            obj.Cis = eph.Cis;
            obj.deltaN = eph.Delta_n;
            obj.M0 = eph.M0;
            obj.eccentricity = eph.Eccentricity;
            obj.sqrtA = eph.sqrtA;
            obj.toe = eph.Toe;
            obj.omega0 = eph.OMEGA0;
            obj.i0 = eph.i0;
            obj.omegaDot = eph.OMEGA_DOT;
            obj.iDot = eph.IDOT;
            obj.ID = eph.SatelliteID;
        end

        function calcSVPositions(obj,transmitTime,transitTime)

            dt = obj.check_t(transmitTime - obj.transmissionTime);

            satClkCorr = (obj.svClockDriftRate.*dt + obj.svClockDrift).*dt + obj.svClockBias - obj.TGD;

            time = transmitTime - satClkCorr;

            a = obj.sqrtA.^2;

            tk = obj.check_t(time - obj.toe);

            n0 = sqrt(obj.GM/a.^3);

            n = n0 + obj.deltaN;

            M = obj.M0 + n*tk;

            M = rem(M + 2*obj.gpsPi, 2*obj.gpsPi);

            E = M;

            for ii = 1:10
                E = M + obj.eccentricity.*sin(E);
            end

            E = rem(E + 2*obj.gpsPi, 2*obj.gpsPi);

            dtr = obj.F.*obj.eccentricity.*obj.sqrtA.*sin(E);

            nu = atan2(sqrt(1 - obj.eccentricity.^2).*sin(E), cos(E) - obj.eccentricity);

            phi = nu + obj.omega0;

            phi = rem(phi, 2*obj.gpsPi);

            u = phi + obj.Cuc.*cos(2*phi) + obj.Cus.*sin(2*phi);

            r = a.*(1 - obj.eccentricity.*cos(E)) + obj.Crc.*cos(2*phi) + obj.Crs.*sin(2*phi);

            incl = obj.i0 + obj.iDot.*tk + obj.Cic.*cos(2*phi) + obj.Cis.*sin(2*phi);

            Omega = obj.omega0 + (obj.omegaDot - obj.Omegae_dot).*tk - obj.Omegae_dot.*obj.toe - obj.Omegae_dot.*transitTime;
        
            Omega = rem(Omega + 2*obj.gpsPi, 2*obj.gpsPi);

            X = r.*cos(u);
            Y = r.*sin(u);

            x = X.*cos(Omega) - Y.*cos(incl).*sin(Omega);
            y = X.*sin(Omega) + Y.*cos(incl).*cos(Omega);
            z = Y.*sin(incl);

            Edot = (n0' + obj.deltaN)./(1-obj.eccentricity.*cos(E));

            phidot = (sqrt(1 - obj.eccentricity.^2)./(1 - obj.eccentricity.*cos(E))).*Edot;

            udot = (1 + 2.*obj.Cus.*cos(2*phi) - 2*obj.Cuc.*sin(2*phi)).*phidot;

            rdot = 2*(obj.Crs.*cos(2*phi)-obj.Crc.*sin(2*phi)).*phidot + a.*obj.eccentricity.*sin(E).*Edot;

            idot = 2*(obj.Cis.*cos(2*phi) - obj.Cic.*sin(2*phi)).*phidot + obj.iDot;

            Xdot = rdot.*cos(u) - r.*sin(u).*udot;

            Ydot = rdot.*sin(u) + r.*cos(u).*udot;

            deltaOmegaDot = obj.omegaDot - obj.Omegae_dot;

            xdot = Xdot.*cos(Omega) - Ydot.*cos(incl).*sin(Omega) + Y.*sin(incl).*sin(Omega).*idot - y.*deltaOmegaDot;

            ydot = Xdot.*sin(Omega) + Ydot.*cos(incl).*cos(Omega) - Y.*sin(incl).*cos(Omega).*idot + x.*deltaOmegaDot;

            zdot = Ydot.*sin(incl) + Y.*cos(incl).*idot;

            obj.svPosX = x;
            obj.svPosY = y;
            obj.svPosZ = z;

            obj.svVelX = xdot;
            obj.svVelY = ydot;
            obj.svVelZ = zdot;

            obj.svClockCorr = ((obj.svClockDriftRate.*dt + obj.svClockDrift).*dt + obj.svClockBias - obj.TGD + dtr)';

        end

        function corrTime = check_t(~,time)

            half_week = 302400;

            corrTime = time;

            if time > half_week
                corrTime = time - 2*half_week;
            elseif time < -half_week
                corrTime = time + 2*half_week;
            end
        end

        function calcUnitVectors(obj,userPos)

            svPos = [obj.svPosX obj.svPosY obj.svPosZ];
            deltaPos = svPos - userPos;
            range = vecnorm(deltaPos,2,2);

            uX = deltaPos(:,1)./range;
            uY = deltaPos(:,2)./range;
            uZ = deltaPos(:,3)./range;

            obj.unitVectors = [uX uY uZ];



        end
    end
end
