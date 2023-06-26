classdef SatellitePositions < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties

    end

    methods
        function obj = SatellitePositions(eph,transmitTime,transitTime)



            obj.calcSVPositions(eph,transmitTime,transitTime);


        end

    end
    methods (Access = private)

        function calcSVPositions(obj,eph,transmitTime,transitTime)

            T_GD=eph.TGD;
            t_oc=eph.TransmissionTime;
            a_f2=eph.SVClockDriftRate;
            a_f1=eph.SVClockDrift;
            a_f0=eph.SVClockBias;
            C_rc=eph.Crc;
            C_rs=eph.Crs;
            C_uc=eph.Cuc;
            C_us=eph.Cus;
            C_ic=eph.Cic;
            C_is=eph.Cis;
            Delta_n=eph.Delta_n;
            M_0=eph.M0;
            e=eph.Eccentricity;
            sqrt_A=eph.sqrtA; % Note this code assumes the data is A!!
            t_oe=eph.Toe;
            Omega_0=eph.OMEGA0;
            i_0=eph.i0;
            omega=eph.omega;
            dot_Omega=eph.OMEGA_DOT;
            Idot=eph.IDOT;
            % GPS constatns

            gpsPi = 3.1415926535898;  % Pi used in the GPS coordinate system

            %--- Constants for satellite position calculation -------------------------
            Omegae_dot = 7.2921151467e-5;  % Earth rotation rate, [rad/s]
            GM = 3.986005e14;      % Earth's universal [m^3/s^2]
            F = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]


            %--- Find time difference ---------------------------------------------
            dt = obj.check_t(transmitTime - t_oc);

            %--- Calculate clock correction ---------------------------------------
            satClkCorr = (a_f2.*dt + a_f1).*dt + a_f0 - T_GD;

            time = transmitTime - satClkCorr;

            % Find satellite's position ----------------------------------------------

            %Restore semi-major axis (this is not needed if ephemeris is already given
            % as "a" - In that case simply use:  a=semiMajAx);
            a   = sqrt_A.*sqrt_A;

            %Time correction
            tk  = obj.check_t(time - t_oe);

            %Initial mean motion
            n0  = sqrt(GM/a.^3);
            %Mean motion
            n   = n0 + Delta_n;

            %Mean anomaly
            M   = M_0 + n * tk;

            %Reduce mean anomaly to between 0 and 360 deg
            M   = rem(M + 2*gpsPi, 2*gpsPi);

            %Initial guess of eccentric anomaly
            E   = M;

            %--- Iteratively compute eccentric anomaly ----------------------------
            for ii = 1:10
                E_old = E;
                E = M + e.*sin(E);
                dE = rem(E - E_old, 2*gpsPi);
            end

            %Reduce eccentric anomaly to between 0 and 360 deg
            E = rem(E + 2*gpsPi, 2*gpsPi);

            %Compute relativistic correction term
            dtr = F.*e.*sqrt_A.*sin(E);

            %Calculate the true anomaly
            nu = atan2(sqrt(1 - e.^2).*sin(E), cos(E)-e);

            %Compute angle phi
            phi = nu + omega;

            %Reduce phi to between 0 and 360 deg
            phi = rem(phi, 2*gpsPi);

            %Correct argument of latitude
            u = phi + C_uc.*cos(2*phi) + C_us.*sin(2*phi);
            %Correct radius
            r = a.*(1 - e.*cos(E)) + C_rc.*cos(2*phi) + C_rs.*sin(2*phi);
            %Correct inclination
            i = i_0 + Idot.*tk + C_ic.*cos(2*phi) + C_is.*sin(2*phi);

            %Compute the angle between the ascending node and the Greenwich meridian
            Omega = Omega_0 + (dot_Omega - Omegae_dot).*tk - Omegae_dot.*t_oe - Omegae_dot.*transitTime;
            %Reduce to between 0 and 360 deg
            Omega = rem(Omega + 2*gpsPi, 2*gpsPi);

            X = r.*cos(u);
            Y = r.*sin(u);

            %--- Compute satellite coordinates ------------------------------------
            x = X.*cos(Omega) - Y.*cos(i).*sin(Omega);
            y = X.*sin(Omega) + Y.*cos(i).*cos(Omega);
            z = Y.*sin(i);

            %% SV velocity calculations
            Edot = (n0' + Delta_n)./(1-e.*cos(E));
            phidot = (sqrt(1 - e.^2)./(1 - e.*cos(E))).*Edot;
            udot = (1 + 2.*C_us.*cos(2*phi) - 2*C_uc.*sin(2*phi)).*phidot;
            rdot = 2*(C_rs.*cos(2*phi)-C_rc.*sin(2*phi)).*phidot + a.*e.*sin(E).*Edot;
            idot = 2*(C_is.*cos(2*phi) - C_ic.*sin(2*phi)).*phidot + Idot;
            Xdot = rdot.*cos(u) - r.*sin(u).*udot;
            Ydot = rdot.*sin(u) + r.*cos(u).*udot;
            Omegadot = dot_Omega - Omegae_dot;
            xdot = Xdot.*cos(Omega) - Ydot.*cos(i).*sin(Omega) + Y.*sin(i).*sin(Omega).*idot - y.*Omegadot;
            ydot = Xdot.*sin(Omega) + Ydot.*cos(i).*cos(Omega) - Y.*sin(i).*cos(Omega).*idot + x.*Omegadot;
            zdot = Ydot.*sin(i) + Y.*cos(i).*idot;

            %%


            %--- Compute satellite coordinates ------------------------------------
            satPositions(1,:) = x;
            satPositions(2,:) = y;
            satPositions(3,:) = z;

            satVelocities(1,:) = xdot;
            satVelocities(2,:) = ydot;
            satVelocities(3,:) = zdot;


            % Include relativistic correction in clock correction -----------------
            satClkCorr = ((a_f2.*dt + a_f1).*dt + a_f0 - T_GD + dtr)';

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
    end
end
