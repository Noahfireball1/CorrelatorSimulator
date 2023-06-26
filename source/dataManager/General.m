classdef General < handle
    %GENERAL Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
        year = [];
        month = [];
        day = [];
        time = [];
        verbose = [];
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
            switch obj.extension
                case ".rnx.gz"
                    eph = rinexread(rinexFilePath);
                    obj.ephemeris = eph.GPS;

                case ".Z"

            end

        end

        function rinexFile = loadRinexFile(obj,dir)

            yearSuffix = num2str(obj.year);

            baseURL = "https://cddis.nasa.gov/archive/gnss/data/daily/";
            webURL = sprintf('%s/%s/%sn/',num2str(obj.year),string(obj.julianDay),yearSuffix(3:end));
            fileURL = sprintf('AMC400USA_R_%s%s0000_01D_GN.rnx.gz',num2str(obj.year),string(obj.julianDay));
            url = append(baseURL, webURL,fileURL);
            outputFilePath = fullfile(dir.data,fileURL);
            [path,unzippedFile,~] = fileparts(outputFilePath);

            if ~exist(append(path,filesep,unzippedFile),'file')

                obj.consoleText(1)

                obj.downloadRinex(url,dir,fileURL)

                obj.uncompressRinex(dir,fileURL)

                obj.consoleText(2)

            end

            rinexFile = append(path,filesep,unzippedFile);

        end

        function downloadRinex(obj,url,dir,fileURL)

            obj.checkPassword(dir);

            system(sprintf('curl -s -c %scookies.txt --ciphers DEFAULT@SECLEVEL=1 --netrc-file %s -L -O %s',dir.data,dir.netRCFilePath,url));
            
            system(sprintf('mv %s %s',fileURL,dir.data));

            system(sprintf('del %scookies.txt',dir.data));
            
        end

        function uncompressRinex(obj,dir,fileURL)

            switch obj.extension
                case ".rnx.gz"
                    gunzip(append(dir.data,fileURL));
                    system(sprintf('del %s%s',dir.data,fileURL));
                case ".gz"
                    uncompress(append(dir.data,fileURL));
            end

        end

        function checkPassword(~,dir)

            if isempty(dir.login)

                fprintf('[correlator-sim] https://cddis.nasa.gov username and password not found.\nPlease enter below:\n')
                username = input('Username:','s');
                password = input('Password:','s');
                dir().createCDDISLoginFile(username,password);
                dir().createNetRCFile(username,password);
            elseif ~exist(append(dir.data,'.netrc'),"file")

                username = dir.login.username;
                password = dir.login.password;

                dir().createNetRCFile(username,password);
            end
        end

        function consoleText(obj,num)

            if obj.verbose
                switch num
                    case 1
                        fprintf('[correlator-sim] Rinex file for specified date not found. Downloading from https://cddis.nasa.gov\n')
                    case 2
                        fprintf('[correlator-sim] Rinex File successfully downloaded and unzipped! \n')
                    case 3

                    otherwise
                end
            end
        end
    end

    methods

        function julianDay = get.julianDay(obj)
            julianDay = day(datetime(obj.year,obj.month,obj.day),"dayofyear");
        end

        function extension = get.extension(obj)
            if obj.year > 2020
                extension = ".rnx.gz";
            else
                extension = ".gz";
            end
        end

    end

end