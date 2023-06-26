classdef Trajectory
    %TRAJECTORY defines trajectory properties that the user defines

    properties
        motion = [];
        latitude = [];
        longitude = [];
        altitude = [];
        position = [];
        velocity = [];
        clockBias = [];
        clockDrift = [];
        ionosphereDelay = [];
        troposphereDelay = [];
    end

    methods
        function obj = Trajectory(config)

            obj.motion = config.motion;
            obj.latitude = config.latitude;
            obj.longitude = config.longitude;
            obj.altitude = config.altitude;
            obj.position = lla2ecef([obj.latitude obj.longitude obj.altitude]);
            obj.velocity = str2num(config.velocity{1});
            obj.clockBias = config.clockBias;
            obj.clockDrift = config.clockDrift;
            obj.ionosphereDelay = config.ionosphereDelay;
            obj.troposphereDelay = config.troposphereDelay;

        end
    end
end