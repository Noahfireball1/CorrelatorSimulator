classdef PropagateCorrelators < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (Constant)
        C = 299792458;
    end

    methods
        function obj = PropagateCorrelators(sim,sampleCount)

            for sv = 1:sim.numChannels

                channel = sprintf('channel%i',sv);

                % Update Reference Recieve Time
                sim.reference.(channel).receiveTime = sim.reference.initialReceiveTime + sampleCount/obj.C;

                % Update Estimated Receive Time
                sim.estimate.(channel).receiveTime = sim.estimate.(channel).initialReceiveTime + sampleCount/obj.C;

                previousSVPos = [sim.satellitePositions.svPosX(sv)...
                    sim.satellitePositions.svPosY(sv)...
                    sim.satellitePositions.svPosZ(sv)];

                % Calculate Previous Time Step's Range
                previousRange = norm(previousSVPos - sim.traj.position);

                % Calculate Transit Time of Signal
                transitTime = previousRange/obj.C;

                % Calculate Satellite Position
                for i = 1:10

                    % Calculate Satellite Positions
                    transmitTime = sim.reference.(channel).receiveTime - transitTime;
                    [x,y,z,u,v,w,clk] = obj.calcSVPositions(obj,transmitTime,transitTime);

                    % Recalculate Transit Time
                    svPos = [x,y,z];
                    transitTime = vecnorm((svPos - position),2,2)/obj.C;

                end
                sim.satellitePositions.svPosX(sv) = x;
                sim.satellitePositions.svPosY(sv) = y;
                sim.satellitePositions.svPosZ(sv) = z;
                sim.satellitePositions.svVelX(sv) = u;
                sim.satellitePositions.svVelY(sv) = v;
                sim.satellitePositions.svVelZ(sv) = w;
                sim.satellitePositions.svClockCorr(sv) = clk;

                % Calculate Reference Unit Vectors
                [ux,uy,uz] = obj.calcUnitVectors(obj,svX,svY,svZ,userPos);

                % Update Reference Parameters

                % Update Line of Sight Velocity

                % Update Carrier Frequency

                % Update Doppler

                % Update Code Frequency

                % Update Reference Range

                % Update Reference Transmit Time

                % Previous Time Step's Code Phase

                % Update Reference Code Phase

                % Update Carrier Phase Advance

                % Update Carrier Phase

                % Calculate Estimated Time Difference

                % Calculate Estimated Code Phase

                % Update Estimated Transmit Time

                % Truncate Code Phase

                % Calculate Estimated Carrier Phase Advance

                % Calculate Estimated Carrier Phase

                % Update Accumulated Doppler

                % Calculate Code Phase Error

                % Accumulate Code Phase Error

                % Accumulate Carrier Phase Error

                % Accumulate Carrier Frequency Error

                % Increment Samples Counter

                % Closing Tracking Loop When Conditions are Met


                stop = 1;
            end
            stop = 1;
        end
    end

    methods (Access = private)
        function [x,y,z,u,v,w,clk] = calcSVPositions(obj,transmitTime,transitTime)

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

        function [ux,uy,uz] = calcUnitVectors(~,svX,svY,svZ,userPos)

            svPos = [svX,svY,svZ];
            deltaPos = svPos - userPos;
            range = vecnorm(deltaPos,2,2);

            ux = deltaPos(:,1)./range;
            uy = deltaPos(:,2)./range;
            uz = deltaPos(:,3)./range;
        end
    end
end