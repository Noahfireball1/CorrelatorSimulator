classdef UserDirectories < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        project = [];
        config = [];
        data = [];
        output = [];
        source = [];
        login = []
        netRCFilePath = [];
    end

    methods
        function obj = UserDirectories()

            obj.project = fileparts(which("run_simulation.m"));
            obj.config = append(obj.project,filesep,'config',filesep);
            obj.data = append(obj.project,filesep,'data',filesep);
            obj.output = append(obj.project,filesep,'output',filesep);
            obj.source = append(obj.project,filesep,'source',filesep);
            obj.readLoginTextFile();

            if exist(append(obj.data,'.netrc'),"file")
                obj.netRCFilePath = append(obj.data,'.netrc');
            end


        end

        function readLoginTextFile(obj)

            if exist(append(obj.config,'CDDISLogin.txt'),"file")
                textFilePath = append(obj.config,'CDDISLogin.txt');
                textFile = fileread(textFilePath);
                splitLogin = strsplit(textFile,{':','\n'});

                obj.login.username = splitLogin{2}(1:end-1);
                obj.login.password = splitLogin{4};
            end


        end

        function createCDDISLoginFile(obj,username,password)

            fid = fopen(append(obj.config,'CDDISLogin.txt'),"wt");
            fprintf(fid,'username:%s\npassword:%s\n',username,password);
            fclose(fid);

        end

        function createNetRCFile(obj,username,password)

            fid = fopen(append(obj.data,'.netrc'),"wt");
            fprintf(fid,'machine urs.earthdata.nasa.gov login %s password %s',username,password);
            fclose(fid);

            obj.netRCFilePath = append(obj.data,'.netrc');

        end

    end
end