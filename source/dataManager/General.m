classdef General < handle
    %GENERAL Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
        year = [];
        month = [];
        day = [];
        time = [];
        verbose = [];
    end

    properties (Access = protected)
        ephemeris = [];
    end

    properties (Dependent)

        julianDay;
        extension;

    end


    methods (Access = public)
        function obj = General(config,dir)
            obj.year = config.year;
            obj.month = config.month;
            obj.day = config.day;
            obj.time = config.time;
            obj.verbose = config.verbose;

            obj.loadEphemeris(dir);
        end
    end

    methods (Access = private)
        function loadEphemeris(obj,dir)

            % Load or Download Correct Rinex File
            rinexFilePath = obj.loadRinexFile(dir);

            % Parse Ephemeris from File



        end

        function rinexFile = loadRinexFile(obj,dir)

            yearSuffix = num2str(obj.year);

            baseURL = "https://cddis.nasa.gov/archive/gnss/data/daily/";
            fileURL = sprintf('%s/%s/%sn/',num2str(obj.year),string(obj.julianDay),yearSuffix(3:end));
            filePath = sprintf('%s_%s_%sn',num2str(obj.year),string(obj.julianDay),yearSuffix(3:end));
            url = append(baseURL, fileURL);
            output_file_path = fullfile(dir.data,filePath);
            [output_parent, output_name, output_ext] = fileparts(output_file_path);

            if ~exist(output_file_path,'file')

                cddis_request(url, output_file_path)

                switch obj.extension
                    case "gz"
                        gunzip(output_file_path)
                    case "Z"
                        uncompress(output_file_path)
                end

                fprintf("Successfully downloaded and unzipped the following file from cddis.nasa.gov: %s\n", ...
                    strcat(output_name, output_ext))
            end

            unzipped_file_path = fullfile(output_parent, output_name);

        end

        function uncompress(file_name)

            if ispc
                str = strjoin({'7z e', char(file_name)});
            else
                str = strcat('uncompress', '  ', file_name);
            end

            system(string(str));
        end



    end

    methods

        function julianDay = get.julianDay(obj)
            julianDay = day(datetime(obj.year,obj.month,obj.day),"dayofyear");
        end

        function extension = get.extension(obj)
            if obj.year > 2020
                extension = "gz";
            else
                extension = "Z";
            end
        end

    end

end