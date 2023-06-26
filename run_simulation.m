%% Formatting
clc
clear
close all
format shortg

useGUI = 0;

%% Adding Paths Based on User's Directories
projectRoot = fileparts(which(mfilename));
addpath(genpath(projectRoot))
configDir = append(projectRoot,filesep,'config',filesep);
dataDir = append(projectRoot,filesep,'data',filesep);
outputDir = append(projectRoot,filesep,'output',filesep);

%% Select a Configuration File
if useGUI
    inputFile = uigetfile('*.yaml','Select Input File',configDir);
    inputFilePath = append(configDir,inputFile);
else
    inputFilePath = append(configDir,'example.yaml');
end

%% Initializing Simulation
settings = DataManager(inputFilePath);
sim = CorrelatorSim(settings);

%% Calculating Initial Satellite Positions
sim.calcSVPos;
%% Initialize Reference Data
sim.initializeReference;
%% initialize Estimated Data
sim.initializeEstimate;
%% Initialize Errors
% error initialized to zero to reflect perfect tracking data at start of
% simulation
initial_position_error = [0.;0.;0.];
initial_range_error = norm(initial_position_error);
initial_velocity_error = [0.0;0.0;0.0];

%% Initialize Scalar Estimates
% clock drift is initialize to zero (?)
% remaining states match true states
est.xyz = true.xyz - initial_position_error;
[est.lat est.lon est.alt] = wgsxyz2lla(est.xyz);
est.vel.xyz = true.vel.xyz + initial_velocity_error;
est.clock.bias = 0;
est.clock.drift = 0;
est.state.vector = [est.xyz(1);est.vel.xyz(1);...
    est.xyz(2);est.vel.xyz(2);...
    est.xyz(3);est.vel.xyz(3);...
    est.clock.bias;est.clock.drift];
est.state.posIdx = [1;3;5];
est.state.velIdx = [2;4;6];
est.state.clkIdx = 7;
est.state.dftIdx = 8;
est.state.covariance = diag([10000*ones(1,6),10000,10000]);

% Setup Estimated Channel Data
for kk = 1:length(true.prns)

    % prn number be tracked
    est.channel(true.prns(kk)).prn = true.prns(kk);

    % carrier to noise ratio in Hz of received signal
    est.channel(true.prns(kk)).cnoL1 = 10^(45/10);

    % line of sight velocity and resulting carrier frequency
    velData = calculate_frequency(unitVector(:,true.prns(kk)),est.vel.xyz,svvel(:,true.prns(kk)),true.clock.drift,constants);

    % phase error may be calculated as a function of initial position error
    phaseErr = calculate_phase_error(unitVector(:,true.prns(kk)),est.xyz,true.xyz,0,constants);

    % range error as a function of initial position errror
    initial_range_error = phaseErr/constants.chipRate*constants.c;

    % save range error for future analysis
    save_init_range_error(kk) = initial_range_error;

    % record line of sight velocity
    est.channel(true.prns(kk)).losVel = velData(1);

    % record initial carrier frequency
    est.channel(true.prns(kk)).carrFreq = velData(2);

    % record initial doppler frequency
    est.channel(true.prns(kk)).doppler = velData(2) - constants.IF;

    % initial code frequnecy
    est.channel(true.prns(kk)).codeFreq = est.channel(true.prns(kk)).doppler * constants.chipRate/constants.L1 + constants.chipRate;

    % initial code phase error as a fuction of initial position errror
    initial_code_phase_error = initial_range_error*constants.chipRate/constants.c;

    % initial code phase value
    est.channel(true.prns(kk)).codePhase = true.channel(true.prns(kk)).codePhase - phaseErr;

    % initial carrier phase error as a function of initial position error
    initial_carrier_phase_error = mod((initial_range_error*true.channel(true.prns(kk)).carrFreq/constants.c*2.*pi),2*pi);
    if initial_carrier_phase_error > pi
        initial_carrier_phase_error = initial_carrier_phase_error - 2*pi;
    end

    % initial carrier phase value
    est.channel(true.prns(kk)).phase = true.channel(true.prns(kk)).phase - initial_carrier_phase_error;

    % initial tracking loop error signals
    est.channel(true.prns(kk)).codeError = 0;
    est.channel(true.prns(kk)).oldCarrierError = 0;
    est.channel(true.prns(kk)).oldCodeError = 0;
    est.channel(true.prns(kk)).carrierError = 0;
    est.channel(true.prns(kk)).oldFreqError = 0;
    est.channel(true.prns(kk)).oldCodeFreq = est.channel(true.prns(kk)).codeFreq;
    est.channel(true.prns(kk)).oldCarrerFreq = est.channel(true.prns(kk)).carrFreq;

    % tracking loop as kalman filter
    est.kal(true.prns(kk)).A = [1, constants.pdiTime, constants.pdiTime^2/2; 0, 1, constants.pdiTime; 0, 0, 1];
    est.kal(true.prns(kk)).C = [1, constants.pdiTime/2, constants.pdiTime^2/6];
    est.kal(true.prns(kk)).B = [-constants.pdiTime;0;0];
    est.kal(true.prns(kk)).X = [0;est.channel(true.prns(kk)).carrFreq;0];
    est.kal(true.prns(kk)).L = [0.541080410673989;7.83430131462034;56.6176707167131];

    % initial receive time
    est.channel(true.prns(kk)).receiveTime = true.channel(true.prns(kk)).receiveTime;

    % initial transmit time in whole CA code periods (i.e. mod 1e-3)
    est.channel(true.prns(kk)).transmitTime = true.channel(true.prns(kk)).transmitTime - mod(true.channel(true.prns(kk)).transmitTime,1e-3);% + svclk(true.prns(kk));

    % counter for recording tracking parameters
    est.channel(true.prns(kk)).updateCounter = 1;

    % counter for counting CA code periods
    est.channel(true.prns(kk)).codePeriodCounter = 0;

    % record initial receive time
    est.channel(true.prns(kk)).initial_receive_time = est.channel(true.prns(kk)).receiveTime;

    % record initial code phase
    est.channel(true.prns(kk)).initial_code_phase = est.channel(true.prns(kk)).codePhase;

    % record initial transmit time
    est.channel(true.prns(kk)).initial_transmit_time = est.channel(true.prns(kk)).transmitTime;

    % record initial carrier phase measurement
    est.channel(true.prns(kk)).accDop = - est.channel(true.prns(kk)).phase;

    % record initial carrier phase measurement for delta carrier phase
    % range rate determinatino
    est.channel(true.prns(kk)).oldCarr = est.channel(true.prns(kk)).accDop;

    % flag for triggering loop closure
    est.channel(true.prns(kk)).newCodePeriodStarted = 0;

end

%% least square position estimation from carrier phase
lsPos = est.state.vector(est.state.posIdx);
lsVel = est.state.vector(est.state.velIdx);
lsClk = est.state.vector(est.state.clkIdx);
lsDft = est.state.vector(est.state.dftIdx);
ls.state.vector = est.state.vector;
ls.state.covariance = est.state.covariance;

%% Tracking Loop Parameters
% third order PLL
carrierBW = 10;
carrierBW = carrierBW/0.7845;
codeA = 1.1;
codeB = 2.4;
track.third.kCarrier1 = carrierBW^3*(constants.pdiTime^2) + codeA*carrierBW^2*constants.pdiTime + codeB*carrierBW;
track.third.kCarrier2 = -(codeA*carrierBW^2*(constants.pdiTime) + 2*codeB*carrierBW);
track.third.kCarrier3 = codeB*carrierBW;

% second order PLL
% carrierBW = 7.;
carrierBW2 = 2;
carrierWN2 = carrierBW/0.53;
carrierA = 1.414;
track.second.kCarrier1 = carrierWN2^2*(constants.pdiTime) + carrierA*carrierWN2;
track.second.kCarrier2 = -carrierA*carrierWN2;

% third order DLL
codeBW = 2;
codeWN = codeBW/0.7845;
codeA = 1.1;
codeB = 2.4;
track.third.kCode1 = codeWN^3*(constants.pdiTime^2) + codeA*codeWN^2*constants.pdiTime + codeB*codeWN;
track.third.kCode2 = -(codeA*codeWN^2*(constants.pdiTime) + 2*codeB*codeWN);
track.third.kCode3 = codeB*codeWN;

% second order DLL
codeBW = 2.;
codeWN = carrierBW/0.53;
codeA = 1.414;
track.second.kCode1 = codeWN^2*(constants.pdiTime) + codeA*codeWN;
track.second.kCode2 = -codeA*codeWN;

% second order FLL
fllBW = 3.5;
fllWN = carrierBW/0.53;
carrierA = 1.414;
track.second.kFreq1 = fllWN^2*(constants.pdiTime)^2 + carrierA*fllWN*constants.pdiTime;
track.second.kFreq2 = -carrierA*fllWN*constants.pdiTime;

%% Saved Data Structs
% data stucture for storing channel data
scalar.carrierError = NaN(true.numChannels,channel_data_length);
scalar.avgCarrError = NaN(true.numChannels,channel_data_length);
scalar.codeError = NaN(true.numChannels,channel_data_length);
scalar.carrFreqError = NaN(true.numChannels,channel_data_length);
scalar.codeFreqError = NaN(true.numChannels,channel_data_length);
scalar.carrDisc = NaN(true.numChannels,channel_data_length);
scalar.codeDisc = NaN(true.numChannels,channel_data_length);
scalar.IP = NaN(true.numChannels,channel_data_length);
scalar.QP = NaN(true.numChannels,channel_data_length);
scalar.IE = NaN(true.numChannels,channel_data_length);
scalar.QE = NaN(true.numChannels,channel_data_length);
scalar.IL = NaN(true.numChannels,channel_data_length);
scalar.QL = NaN(true.numChannels,channel_data_length);
scalar.receiveTime = NaN(true.numChannels,channel_data_length);
scalar.transmitTime = NaN(true.numChannels,channel_data_length);
scalar.rangeRes = NaN(true.numChannels,channel_data_length);

% data structure for storing estimated navigation solution
nav.pos = NaN(3,nav_data_length);
nav.vel = NaN(3,nav_data_length);
nav.clkBias = NaN(1,nav_data_length);
nav.clkDrft = NaN(1,nav_data_length);
nav.rangeRes = NaN(true.numChannels,nav_data_length);
nav.rateRes = NaN(true.numChannels,nav_data_length);
nav.receiveTime = NaN(1,nav_data_length);
nav.cnoL1 = NaN(true.numChannels,nav_data_length);
nav.lsPos = NaN(3,nav_data_length);
nav.lsVel = NaN(3,nav_data_length);
nav.lsClkBias = NaN(1,nav_data_length);
nav.lsClkDrft = NaN(1,nav_data_length);
nav.covariance = NaN(8,nav_data_length);
nav.rangeResVar = NaN(true.numChannels,nav_data_length);
nav.rateResVar = NaN(true.numChannels,nav_data_length);


% data structure for storing all true parameters channel and nav
reference.carrPhase = NaN(true.numChannels,channel_data_length);
reference.codePhase = NaN(true.numChannels,channel_data_length);
reference.carrFreq = NaN(true.numChannels,channel_data_length);
reference.codeFreq = NaN(true.numChannels,channel_data_length);
reference.receiveTime = NaN(true.numChannels,channel_data_length);
reference.transmitTime = NaN(true.numChannels,channel_data_length);

% true nav data
reference.el = NaN(true.numChannels,nav_data_length);
reference.az = NaN(true.numChannels,nav_data_length);
reference.rangeError = NaN(true.numChannels,nav_data_length);
reference.posError = NaN(3,nav_data_length);
reference.pos = NaN(3,nav_data_length);
reference.vel = NaN(3,nav_data_length);
reference.clkBias = NaN(1,nav_data_length);
reference.clkDrft = NaN(1,nav_data_length);

% vector for intermediate calculations
codePhaseErrAccum = zeros(true.numChannels,1);
phaseErrAccum = zeros(true.numChannels,1);
freqErrAccum = zeros(true.numChannels,1);
integratedSamplesCount = zeros(true.numChannels,1);
avgAmp = zeros(true.numChannels,1);
avgNoiseVar = zeros(true.numChannels,1);
pastMiddle = zeros(true.numChannels,1);
trajectoryUpdate = ones(true.numChannels,1);
measTrueRange = zeros(true.numChannels,1);
lsPhase = zeros(true.numChannels,1);
lsDeltaPhase = zeros(true.numChannels,1);

%% Main Loop Parameters
eof = 0;
count = 1;
firstLoop = 1;
sampleCount = 0;
navUpdateCounter = 1;

% main loop
while ~eof

    % processing next sample (sampleCount initialize to zero)
    sampleCount = sampleCount + 1;

    if ~mod(sampleCount/sample_run_length*100,10)

        % print percentage to screen
        sprintf('Percent Complete = %f',sampleCount/sample_run_length*100)

    end

    % kalman time update
    % calculate state transition matrix and process noise matrix
    [phi, Q] = get_state_transition_process_noise(1/constants.S);

    % propagate states to middle of integration period
    true.state.vector = phi*true.state.vector;

    % update clock parameters
    true.clock.bias = true.state.vector(est.state.clkIdx);
    true.clock.drift = true.state.vector(est.state.dftIdx);

    % update true x y z position
    true.xyz = true.state.vector(est.state.posIdx);

    % update true x y z vel
    true.vel.xyz = true.state.vector(est.state.velIdx);

    % Calculate Correlator Outputs
    for kk = 1:length(true.prns)

        % update true receive time
        true.channel(true.prns(kk)).receiveTime = true.initialReceiveTime + sampleCount/constants.S;

        % udpate estimate receive time
        est.channel(true.prns(kk)).receiveTime = est.channel(true.prns(kk)).initial_receive_time + sampleCount/constants.S;

        % update true signal parameters
        % calculate true transit time of signal
        transitTime = norm(init.svpos(:,true.prns(kk))-true.xyz)/constants.c;

        % old range
        oldRange = norm(init.svpos(:,true.prns(kk))-true.xyz);

        % calculate satellite position at end of integration period
        satPosChange = 1;
        oldsvclk = init.svclk(true.prns(kk));
        while satPosChange > 1e-3
            [svpos(:,true.prns(kk)),svvel(:,true.prns(kk)),svclk(true.prns(kk))] = ...
                calc_sv_pos3(ephemerides2(:,true.prns(kk)),...
                true.channel(true.prns(kk)).receiveTime - transitTime,...
                transitTime);
            transitTime = norm(svpos(:,true.prns(kk))-true.xyz)/constants.c;

            satPosChange = norm(svpos(:,true.prns(kk)) - init.svpos(:,true.prns(kk)));
            init.svpos(:,true.prns(kk)) = svpos(:,true.prns(kk));
            init.svvel(:,true.prns(kk)) = svvel(:,true.prns(kk));
            init.svclk(true.prns(kk)) = svclk(true.prns(kk));

        end

        % calculate true unit vectors at end of integration period
        unitVector = get_unit_vectors(true.xyz,svpos);

        % update true parameters
        velData = calculate_frequency(unitVector(:,true.prns(kk)),true.vel.xyz,svvel(:,true.prns(kk)),true.clock.drift,constants);

        % record los velocity at end of integration period
        true.channel(true.prns(kk)).losVel = velData(1);

        % record carrier frequency at end of integration period
        true.channel(true.prns(kk)).carrFreq = velData(2);

        % record doppler at end of integration period
        true.channel(true.prns(kk)).doppler = velData(2) - constants.IF;

        % record code frequency at end of integration period
        true.channel(true.prns(kk)).codeFreq = true.channel(true.prns(kk)).doppler * constants.chipRate/constants.L1 + constants.chipRate;

        % calculate true range at new sample time
        trueRange = norm(svpos(:,true.prns(kk)) - true.xyz);

        % update true transmit time
        true.channel(true.prns(kk)).transmitTime = true.channel(true.prns(kk)).receiveTime - (trueRange + true.ionoDelay + true.tropDelay + true.clock.bias - svclk(true.prns(kk))*constants.c)/constants.c;

        % grab old code phase
        oldCodePhase = true.channel(true.prns(kk)).codePhase;

        % update true received code phase
        true.channel(true.prns(kk)).codePhase = mod(true.channel(true.prns(kk)).transmitTime,1e-3)*constants.chipRate;

        % calculate carrier phase advance
        truePhaseIncrement = constants.IF/constants.S - (trueRange - oldRange)/constants.wavelength;

        % calculate new carrier phase
        true.channel(true.prns(kk)).phase = constants.IF*sampleCount/constants.S - (trueRange - true.ionoDelay + true.tropDelay + true.clock.bias - svclk(true.prns(kk))*constants.c)/constants.wavelength;

        % calculate estimate dt
        estCodePhaseIncrement = est.channel(true.prns(kk)).codeFreq/constants.S;

        % calculate new estimated code phase
        est.channel(true.prns(kk)).codePhase = est.channel(true.prns(kk)).codePhase + estCodePhaseIncrement;

        % update transmit time
        if (est.channel(true.prns(kk)).codePhase > 1023)
            est.channel(true.prns(kk)).newCodePeriodStarted = 1;
            est.channel(true.prns(kk)).codePeriodCounter = est.channel(true.prns(kk)).codePeriodCounter + 1;
            est.channel(true.prns(kk)).transmitTime = est.channel(true.prns(kk)).initial_transmit_time + est.channel(true.prns(kk)).codePeriodCounter*1e-3;
        end

        % truncate code phase
        est.channel(true.prns(kk)).codePhase = mod(est.channel(true.prns(kk)).codePhase,1023);

        % calculate carrier phase advance
        estPhaseIncrement = est.channel(true.prns(kk)).carrFreq/constants.S;

        % calculate new estimated carrier phase
        est.channel(true.prns(kk)).phase = est.channel(true.prns(kk)).phase + estPhaseIncrement;

        % update accumulated doppler
        est.channel(true.prns(kk)).accDop = est.channel(true.prns(kk)).accDop + (constants.IF/constants.S - estPhaseIncrement);

        % get code phase error
        tcodePhaseErr = true.channel(true.prns(kk)).codePhase - est.channel(true.prns(kk)).codePhase;

        if(tcodePhaseErr > 1000)
            tcodePhaseErr = tcodePhaseErr - 1023;
        elseif(tcodePhaseErr < -1000)
            tcodePhaseErr = tcodePhaseErr + 1023;
        end

        % accumulate code phase error
        codePhaseErrAccum(kk) = codePhaseErrAccum(kk) + tcodePhaseErr;

        % check phase error for wrap issue
        tcarrierPhaseErr = true.channel(true.prns(kk)).phase - est.channel(true.prns(kk)).phase;

        % accumulate carrier phase error
        phaseErrAccum(kk) = phaseErrAccum(kk) + tcarrierPhaseErr;

        % accumulate carrier frequency error
        freqErrAccum(kk) = freqErrAccum(kk) + (true.channel(true.prns(kk)).carrFreq - est.channel(true.prns(kk)).carrFreq);

        % increment counter
        integratedSamplesCount(kk) = integratedSamplesCount(kk) + 1;

        % close tracking loop at pdiTimeCodePeriods number of prn code
        % cycles
        if ~mod(est.channel(true.prns(kk)).codePeriodCounter,constants.pdiTimeCodePeriods) && (est.channel(true.prns(kk)).newCodePeriodStarted)

            % used to ensure that statement isn't true on first integration
            % period after loop closure
            est.channel(true.prns(kk)).newCodePeriodStarted = 0;

            % calculate average carrier phase error
            carrPhaseErr = 2*pi*phaseErrAccum(kk)/integratedSamplesCount(kk);

            % calculate average carrier frequency error
            carrFreqErr = freqErrAccum(kk)/integratedSamplesCount(kk);

            % get code phase error
            codePhaseErr = codePhaseErrAccum(kk)/integratedSamplesCount(kk);

            % calculate noise vector for correlator outputs
            noise = randn(6,1);
            if noNoise
                % zero noise if needed
                noise = zeros(6,1);
            end

            % calculate correlator outputs
            IP = sqrt(2*true.channel(true.prns(kk)).cnoL1*constants.pdiTime)*...
                sinc((carrFreqErr)*constants.pdiTime*2*pi)*...
                (1 - abs(codePhaseErr))*...
                cos(carrPhaseErr) + noise(1);
            QP = sqrt(2*true.channel(true.prns(kk)).cnoL1*constants.pdiTime)*...
                sinc((carrFreqErr)*constants.pdiTime*2*pi)*...
                (1 - abs(codePhaseErr))*...
                sin(carrPhaseErr) + noise(2);
            IE = sqrt(2*true.channel(true.prns(kk)).cnoL1*constants.pdiTime)*...
                sinc((carrFreqErr)*constants.pdiTime*2*pi)*...
                (1 - abs(codePhaseErr - constants.offset))*...
                cos(carrPhaseErr) + noise(3);
            QE = sqrt(2*true.channel(true.prns(kk)).cnoL1*constants.pdiTime)*...
                sinc((carrFreqErr)*constants.pdiTime*2*pi)*...
                (1 - abs(codePhaseErr - constants.offset))*...
                sin(carrPhaseErr) + noise(4);
            IL = sqrt(2*true.channel(true.prns(kk)).cnoL1*constants.pdiTime)*...
                sinc((carrFreqErr)*constants.pdiTime*2*pi)*...
                (1 - abs(codePhaseErr + constants.offset))*...
                cos(carrPhaseErr) + noise(5);
            QL = sqrt(2*true.channel(true.prns(kk)).cnoL1*constants.pdiTime)*...
                sinc((carrFreqErr)*constants.pdiTime*2*pi)*...
                (1 - abs(codePhaseErr + constants.offset))*...
                sin(carrPhaseErr) + noise(6);

            % calculate discriminators
            % scalar carrier phase discriminator -- 2 quadrant arc tangent --
            % units cycles
            carrDisc = atan(QP/IP)/2/pi;

            % scalar code phase discriminator -- normalized early minus late
            % power -- units chips
            codeDisc = 1/4*((IE^2+QE^2) - (IL^2 + QL^2))/((IE^2+QE^2) + (IL^2 + QL^2));

            if (imag(codeDisc))
                stop = 1;
            end

            % maintain estimates of amplitude and noise
            % calculate signal plus noise amplitude
            amp = (IE + IL)^2 + (QE + QL)^2;
            % calcualte noise variance
            noiseVar = var(randn(18,1));
            % running average on amplitude and noise
            if firstLoop
                avgAmp(kk) = amp;
                avgNoiseVar(kk) = noiseVar;
                if(kk==length(true.prns))
                    firstLoop = 0;
                end
            else
                avgAmp(kk) = 9/10*avgAmp(kk) + 1/10*amp;
                avgNoiseVar(kk) = 9/10*avgNoiseVar(kk) + 1/10*noiseVar;
            end

            % reset accumulaters
            phaseErrAccum(kk) = 0;
            freqErrAccum(kk) = 0;
            integratedSamplesCount(kk) = 0;
            codePhaseErrAccum(kk) = 0;

            % record estimate tracking loop data
            scalar.carrierError(kk,est.channel(true.prns(kk)).updateCounter) = true.channel(true.prns(kk)).phase - est.channel(true.prns(kk)).phase;
            scalar.avgCarrError(kk,est.channel(true.prns(kk)).updateCounter) = carrPhaseErr;
            scalar.codeError(kk,est.channel(true.prns(kk)).updateCounter) = codePhaseErr;
            scalar.carrFreqError(kk,est.channel(true.prns(kk)).updateCounter) = true.channel(true.prns(kk)).carrFreq - est.channel(true.prns(kk)).carrFreq;
            scalar.codeFreqError(kk,est.channel(true.prns(kk)).updateCounter) = true.channel(true.prns(kk)).codeFreq - est.channel(true.prns(kk)).codeFreq;
            scalar.carrDisc(kk,est.channel(true.prns(kk)).updateCounter) = carrDisc;
            scalar.codeDisc(kk,est.channel(true.prns(kk)).updateCounter) = codeDisc;
            scalar.IP(kk,est.channel(true.prns(kk)).updateCounter) = IP;
            scalar.QP(kk,est.channel(true.prns(kk)).updateCounter) = QP;
            scalar.IE(kk,est.channel(true.prns(kk)).updateCounter) = IE;
            scalar.QE(kk,est.channel(true.prns(kk)).updateCounter) = QE;
            scalar.IL(kk,est.channel(true.prns(kk)).updateCounter) = IL;
            scalar.QL(kk,est.channel(true.prns(kk)).updateCounter) = QL;
            scalar.receiveTime(kk,est.channel(true.prns(kk)).updateCounter) = est.channel(true.prns(kk)).receiveTime;
            scalar.transmitTime(kk,est.channel(true.prns(kk)).updateCounter) = est.channel(true.prns(kk)).transmitTime + est.channel(true.prns(kk)).codePhase/constants.chipRate;

            % record true tracking parameters
            reference.codePhase(kk,est.channel(true.prns(kk)).updateCounter) = true.channel(true.prns(kk)).codePhase;
            reference.carrPhase(kk,est.channel(true.prns(kk)).updateCounter) = true.channel(true.prns(kk)).phase;
            reference.codeFreq(kk,est.channel(true.prns(kk)).updateCounter) = true.channel(true.prns(kk)).codeFreq;
            reference.carrFreq(kk,est.channel(true.prns(kk)).updateCounter) = true.channel(true.prns(kk)).carrFreq;
            reference.receiveTime(kk,est.channel(true.prns(kk)).updateCounter) = true.channel(true.prns(kk)).receiveTime;
            reference.transmitTime(kk,est.channel(true.prns(kk)).updateCounter) = true.channel(true.prns(kk)).transmitTime;
            est.channel(true.prns(kk)).updateCounter = est.channel(true.prns(kk)).updateCounter + 1;

            % calculate carrier frequency feedback term from velocity
            % feedback
            % calculate current transmit time
            transmitTime = est.channel(true.prns(kk)).transmitTime + est.channel(true.prns(kk)).codePhase/constants.chipRate;

            % get pseudorange
            psrL1 = (est.channel(true.prns(kk)).receiveTime - transmitTime)*constants.c;

            % calculate true transit time of signal
            transitTime = norm(init.svpos(:,true.prns(kk)) - est.state.vector(est.state.posIdx))/constants.c;

            % calculate satellite position at measurment time
            [esvpos,esvvel,esvclk] = ...
                calc_sv_pos3(ephemerides2(:,true.prns(kk)),...
                transmitTime,...
                transitTime);

            % calculate estimated range
            range = norm(esvpos - true.state.vector(est.state.posIdx));

            % calculate unit vector
            uv = (esvpos - true.state.vector(est.state.posIdx))./range;

            % velocity estimate
            feedbackVel = true.state.vector(est.state.velIdx) + feedbackVel_std*randn(3,1);

            % line of sight vel
            velData = calculate_frequency(uv,feedbackVel,esvvel,true.state.vector(est.state.dftIdx),constants);

            % frequency feedback
            feedbackFreq = velData(2);

            if imag(feedbackFreq)
                stop = 1;
            end

            % frequency error
            feedbackFreqErr = feedbackFreq - est.channel(true.prns(kk)).carrFreq;

            % third order loop filters
            % code phase DLL
            %             tmpCodeFreq = est.channel(true.prns(kk)).codeFreq;
            %             est.channel(true.prns(kk)).codeFreq = (2*est.channel(true.prns(kk)).codeFreq - est.channel(true.prns(kk)).oldCodeFreq) + track.third.kCode1*codeDisc + track.third.kCode2*est.channel(true.prns(kk)).codeError + track.third.kCode3*est.channel(true.prns(kk)).oldCodeError;
            est.channel(true.prns(kk)).codeFreq = est.channel(true.prns(kk)).codeFreq + track.second.kCode1*codeDisc + track.second.kCode2*est.channel(true.prns(kk)).codeError;
            est.channel(true.prns(kk)).oldCodeError = est.channel(true.prns(kk)).codeError;
            est.channel(true.prns(kk)).codeError = codeDisc;
            %             est.channel(true.prns(kk)).oldCodeFreq = tmpCodeFreq;

            if ~feedback_flag
                % carrier phase PLL
                tmpCarrFreq = est.channel(true.prns(kk)).carrFreq;
                est.channel(true.prns(kk)).carrFreq = (2*est.channel(true.prns(kk)).carrFreq - est.channel(true.prns(kk)).oldCarrerFreq) + track.third.kCarrier1*carrDisc + track.third.kCarrier2*est.channel(true.prns(kk)).carrierError + track.third.kCarrier3*est.channel(true.prns(kk)).oldCarrierError;
                est.channel(true.prns(kk)).oldCarrierError = est.channel(true.prns(kk)).carrierError;
                est.channel(true.prns(kk)).carrierError = carrDisc;
                est.channel(true.prns(kk)).oldCarrerFreq = tmpCarrFreq;

            else
                % third order pll with 2nd order fll with velocity feedback
                tmpCarrFreq = est.channel(true.prns(kk)).carrFreq;
                est.channel(true.prns(kk)).carrFreq = tmpCarrFreq ...
                    + track.second.kCarrier1*carrDisc ...
                    + track.second.kCarrier2*est.channel(true.prns(kk)).oldCarrierError ...
                    + track.second.kFreq1*feedbackFreqErr ...
                    + track.second.kFreq2*est.channel(true.prns(kk)).oldFreqError;
                est.channel(true.prns(kk)).oldCarrierError = est.channel(true.prns(kk)).carrierError;
                est.channel(true.prns(kk)).carrierError = carrDisc;
                est.channel(true.prns(kk)).oldFreqError = feedbackFreqErr;

            end

        end

    end

    if sampleCount > (constants.navUpdateTimeSamps*navUpdateCounter)

        % kalman time update
        % calculate state transition matrix and process noise matrix
        [phi, Q] = get_state_transition_process_noise(constants.navUpdateTimeSec);

        % propagate states to middle of integration period
        est.state.vector = phi*est.state.vector;

        % propagate covariance matrix
        est.state.covariance = phi*est.state.covariance*phi' + Q;

        idx = 1;
        for kk = 1:length(true.prns)

            % calculate carrier to noise ratio
            if (avgAmp(kk) > 4*avgNoiseVar(kk))
                rawCno = 10*log10(1/2/constants.pdiTime) + 10*log10((avgAmp(kk) - 4*avgNoiseVar(kk))/avgNoiseVar(kk));
            else
                % used to prevent Kalman gain problem
                % in future measurement with CNO this low should be ignored
                rawCno = 15;
            end

            % calculate current transmit time
            transmitTime = est.channel(true.prns(kk)).transmitTime + est.channel(true.prns(kk)).codePhase/constants.chipRate;

            % get pseudorange
            psrL1 = (est.channel(true.prns(kk)).receiveTime - transmitTime)*constants.c;

            % calculate true transit time of signal
            transitTime = norm(init.svpos(:,true.prns(kk)) - est.state.vector(est.state.posIdx))/constants.c;

            % calculate satellite position at measurment time
            [esvpos,esvvel,esvclk] = ...
                calc_sv_pos3(ephemerides2(:,true.prns(kk)),...
                transmitTime,...
                transitTime);

            % calculate estimated range
            range = norm(esvpos - est.state.vector(est.state.posIdx));

            % calculate range error
            rangeResidual(kk) = psrL1 - (norm(esvpos - est.state.vector(est.state.posIdx)) + est.state.vector(est.state.clkIdx) - init.svclk(true.prns(kk))*constants.c);

            % get carrier frequency
            f_meas = est.channel(true.prns(kk)).carrFreq;

            % get doppler frequency
            dopp = f_meas - constants.IF;

            % calculate range rate
            measRate = -dopp*constants.wavelength;

            % get new carrier phase measurment
            newCarr = est.channel(true.prns(kk)).accDop;

            % calculate delta carrier phase over last navigation update
            % period
            deltaCarr = newCarr - est.channel(true.prns(kk)).oldCarr;

            % calcuate range rate from carrier phase
            %             measRate = deltaCarr*constants.wavelength/constants.navUpdateTimeSec;

            % record carrier phas measurement for next update
            est.channel(true.prns(kk)).oldCarr = newCarr;

            % calculate unit vectors
            uv = (esvpos - est.state.vector(est.state.posIdx))./range;

            % calculate estimated range rate
            predRate = -uv'*(est.state.vector(est.state.velIdx) - esvvel) + est.state.vector(est.state.dftIdx);

            % calculate range rate residual
            rateResidual(kk) = measRate - predRate;

            % kalman measurement matrix
            H(idx:idx+1,:) = [-uv(1) 0 -uv(2) 0 -uv(3) 0 1 0;...
                0 -uv(1) 0 -uv(2) 0 -uv(3) 0 1];

            % kalman measurment vector
            z(idx:idx+1) = [rangeResidual(kk);rateResidual(kk)];

            % kalman measurement variance
            [vR,vRR] = get_range_and_range_rate_variance(rawCno);

            % kalman measurement variance matrix
            R(idx,idx) = vR;
            R(idx+1,idx+1) = vRR;
            idx = idx + 2;

            rangeResVar(kk) = vR;
            rateResVar(kk) = vRR;

            % debugging values for true positions
            trueRange(kk) = (norm(svpos(:,true.prns(kk)) - true.xyz) + true.clock.bias + true.ionoDelay + true.tropDelay - svclk(true.prns(kk))*constants.c);
            measTrueRange(kk) = (true.channel(true.prns(kk)).receiveTime - true.channel(true.prns(kk)).transmitTime)*constants.c;

            % save cno value
            nav.cnoL1(kk,navUpdateCounter) = rawCno;

            % calculate satellite position in enufor sky plot
            sned = wgsdiffxyz2diffned(svpos(:,true.prns(kk)) - true.xyz,true.lat,true.lon);
            reference.el(kk,navUpdateCounter) = atan2(-sned(3), sqrt(sned(1)^2 + sned(2)^2))*180/pi;
            reference.az(kk,navUpdateCounter) = atan2(sned(2),sned(1))*180/pi;

        end

        if any(any(isnan(H)))
            break;
        end

        if navUpdateCounter > 5 && cond((H*est.state.covariance*H' + R)) > 1e10
            break
        end

        % calculale kalman gain
        K = est.state.covariance*H'/(H*est.state.covariance*H' + R);

        % update state vector
        est.state.vector = est.state.vector + K*z';

        % update state covariance
        est.state.covariance = (eye(size(est.state.covariance)) - K*H)*est.state.covariance;

        % record estimated state data
        nav.pos(:,navUpdateCounter) = est.state.vector(est.state.posIdx);
        nav.vel(:,navUpdateCounter) = est.state.vector(est.state.velIdx);
        nav.clkBias(navUpdateCounter) = est.state.vector(est.state.clkIdx);
        nav.clkDrft(navUpdateCounter) = est.state.vector(est.state.dftIdx);
        nav.rangeRes(:,navUpdateCounter) = rangeResidual;
        nav.rangeResVar(:,navUpdateCounter) = rangeResVar;
        nav.rateRes(:,navUpdateCounter) = rateResidual;
        nav.rateResVar(:,navUpdateCounter) = rateResVar;
        nav.receiveTime(navUpdateCounter) = est.channel(true.prns(1)).receiveTime;
        nav.variance(:,navUpdateCounter) = diag(est.state.covariance);

        % record true state data
        reference.pos(:,navUpdateCounter) = true.state.vector(est.state.posIdx);
        reference.vel(:,navUpdateCounter) = true.state.vector(est.state.velIdx);
        reference.clkBias(navUpdateCounter) = true.state.vector(est.state.clkIdx);
        reference.clkDrft(navUpdateCounter) = true.state.vector(est.state.dftIdx);
        reference.rangeError(:,navUpdateCounter) = measTrueRange - trueRange';

        % propogate least sqaures state vector for carrier phase
        % positioning
        ls.state.vector = phi*ls.state.vector;
        ls.state.covariance = phi*ls.state.covariance*phi' + Q;
        lsPos = ls.state.vector(est.state.posIdx);
        lsVel = ls.state.vector(est.state.velIdx);
        lsClk = ls.state.vector(est.state.clkIdx);
        lsDft = ls.state.vector(est.state.dftIdx);

        for kk=1:length(true.prns)

            % get sat positions
            % calculate current transmit time
            %             transmitTime = est.channel(true.prns(kk)).transmitTime + est.channel(true.prns(kk)).codePhase/constants.chipRate;
            transmitTime = true.channel(true.prns(kk)).transmitTime;

            % calculate true transit time of signal
            transitTime = norm(init.svpos(:,true.prns(kk)) - lsPos)/constants.c;

            % calculate satellite position at measurment time
            [esvpos,esvvel,esvclk] = ...
                calc_sv_pos3(ephemerides2(:,true.prns(kk)),...
                transmitTime,...
                transitTime);

            % calculate Range Prediction
            lsRangePred = (norm(esvpos - lsPos) + lsClk - esvclk*constants.c);

            % get Range Meas
            lsRangeMeas = est.channel(true.prns(kk)).accDop*constants.wavelength;

            % get measurement vector
            lsZ(kk) = lsRangeMeas - lsRangePred;

            % calculate unit vectors
            uv = (esvpos - lsPos)./lsRangePred;

            % get geometry matrix
            lsH(kk,:) = [-uv(1) 0 -uv(2) 0 -uv(3) 0 1 0];

            % get covariance weighting
            lsR(kk,kk) = R(2*kk,2*kk);

        end

        % calculate kalman gain
        lsK = ls.state.covariance*lsH'/(lsH*ls.state.covariance*lsH' + lsR);

        % calculate least squares correction
        correction = lsK*lsZ';

        % update least square position
        lsPos = lsPos + correction(est.state.posIdx);

        % update velocity
        lsVel = lsVel + correction(est.state.velIdx);

        % update least squares clk bias estimate
        lsClk = lsClk + correction(est.state.clkIdx);

        % update clock drift
        lsDft = lsDft + correction(est.state.dftIdx);

        ls.state.vector(est.state.posIdx) = lsPos;
        ls.state.vector(est.state.velIdx) = lsVel;
        ls.state.vector(est.state.clkIdx) = lsClk;
        ls.state.vector(est.state.dftIdx) = lsDft;

        %  record carrier phase position solution
        nav.lsPos(:,navUpdateCounter) = lsPos;
        nav.lsVel(:,navUpdateCounter) = lsVel;
        nav.lsClkBias(:,navUpdateCounter) = lsClk;
        nav.lsClkDrft(:,navUpdateCounter) = lsDft;

        %         % predict phase at end of next code period
        %         for kk = 1:length(true.prns)
        %
        %             tmpStates = phi*ls.state.vector;
        %             lsPos = ls.state.vector(est.state.posIdx);
        %             lsVel = ls.state.vector(est.state.velIdx);
        %             lsClk = ls.state.vector(est.state.clkIdx);
        %             lsDft = ls.state.vector(est.state.dftIdx);
        %
        %             % get sat positions
        %             % calculate current transmit time
        %             transmitTime = est.channel(true.prns(kk)).transmitTime + est.channel(true.prns(kk)).codePhase/constants.chipRate + 20e-3;
        %
        %             % calculate true transit time of signal
        %             transitTime = norm(init.svpos(:,true.prns(kk)) - lsPos)/constants.c;
        %
        %             % calculate satellite position at end of integration period
        %             satPosChange = 1;
        %             while satPosChange > 1e-3
        %                 [ssvpos,ssvvel,ssvclk] = ...
        %                     calc_sv_pos3(ephemerides2(:,true.prns(kk)),...
        %                     transmitTime,...
        %                     transitTime);
        %                 transitTime = norm(ssvpos-lsPos)/constants.c;
        %
        %                 satPosChange = norm(ssvpos - init.svpos(:,true.prns(kk)));
        %                 init.svpos(:,true.prns(kk)) = ssvpos;
        %                 init.svvel(:,true.prns(kk)) = ssvvel;
        %                 init.svclk(true.prns(kk)) = ssvclk;
        %             end
        %
        %             % calculate Range Prediction
        %             newRangePred = (norm(ssvpos - lsPos) + lsClk - ssvclk*constants.c);
        %
        %             % calculate new carrier phase assumes zero carrier ambiguity
        %             lsPhase(kk) = newRangePred/constants.wavelength;
        %
        %             % phase increment
        %             lsDeltaPhase(kk) = (constants.IF - (lsPhase(kk) - est.channel(true.prns(kk)).accDop)/constants.pdiTime)/constants.S;
        %
        %         end

        % update navUpdateCounter
        navUpdateCounter = navUpdateCounter + 1;

    end

    if sampleCount > sample_run_length
        eof = 1;
    end

end

if ~feedback_flag
    save(sprintf('scalar_threshold_data/cno%d_codebw%d_carrierbw%d',nomCn0,codeBW,carrierBW))
else
    save(sprintf('scalar_threshold_data/cno%d_codebw%d_carrierbw%d_wFeedback_stdx100_%d',nomCn0,codeBW,carrierBW2,feedbackVel_std*100))
end

%% plot scalar data
if plotScalarData
    figure
    hold on
    plot((scalar.receiveTime - scalar.receiveTime(1,1))',scalar.IP','.b')
    plot((scalar.receiveTime - scalar.receiveTime(1,1))',scalar.IE','.g')
    plot((scalar.receiveTime - scalar.receiveTime(1,1))',scalar.IL','.r')
    title('Scalar InPhase E-P-L')
    ylabel('Amplitude')
    xlabel('Time (s)')
    axis([0 10 10 40])

    figure
    hold on
    plot(scalar.QP','.b')
    plot(scalar.QE','.g')
    plot(scalar.QL','.r')
    title('Scalar QuadPhase')

    figure
    plot(scalar.carrierError')
    title('Scalar carrError')

    figure
    plot(scalar.codeError')
    title('Scalar codeError')

    figure
    plot(scalar.carrFreqError')
    title('Scalar carrFreqError')

    figure
    plot(scalar.codeFreqError')
    title('Scalar codeFreqError')

    figure
    plot(scalar.carrDisc')
    title('Scalar carrDisc')

    figure
    plot(scalar.codeDisc')
    title('Scalar codeDisc')

end

if plotNavData

    figure
    plot(reference.pos(1,:) - nav.pos(1,:));
    hold on
    plot(nav.variance(1,:).^0.5,'r');
    plot(-nav.variance(1,:).^0.5,'r');
    title('Position Error X')

    figure
    plot(reference.pos(2,:) - nav.pos(2,:));
    hold on
    plot(nav.variance(3,:).^0.5,'r');
    plot(-nav.variance(3,:).^0.5,'r');
    title('Position Error Y')

    figure
    plot(reference.pos(3,:) - nav.pos(3,:));
    hold on
    plot(nav.variance(5,:).^0.5,'r');
    plot(-nav.variance(5,:).^0.5,'r');
    title('Position Error Z')

    figure
    plot(reference.vel(1,:) - nav.vel(1,:));
    hold on
    plot(nav.variance(2,:).^0.5,'r');
    plot(-nav.variance(2,:).^0.5,'r');
    title('Velocity Error X')

    figure
    plot(reference.vel(2,:) - nav.vel(2,:));
    hold on
    plot(nav.variance(4,:).^0.5,'r');
    plot(-nav.variance(4,:).^0.5,'r');
    title('Velocity Error Y')

    figure
    plot(reference.vel(3,:) - nav.vel(3,:));
    hold on
    plot(nav.variance(6,:).^0.5,'r');
    plot(-nav.variance(6,:).^0.5,'r');
    title('Velocity Error Z')

    figure
    plot(reference.clkBias - nav.clkBias);
    hold on
    plot(nav.variance(7,:).^0.5,'r');
    plot(-nav.variance(7,:).^0.5,'r');
    title('Clock Bias Error')

    figure
    plot(reference.clkDrft - nav.clkDrft);
    hold on
    plot(nav.variance(8,:).^0.5,'r');
    plot(-nav.variance(8,:).^0.5,'r');
    title('Clock Drift Error')

end

if plotVerbose

    figure
    plot(reference.rangeError')
    title('True Psuedorange Time - True Pseudorange Distance')

    figure
    skyPlot(reference.az,reference.el,true.prns,'.')

    figure
    plot((reference.receiveTime - reference.receiveTime(1,1))',(reference.transmitTime - scalar.transmitTime)'*constants.c)
    title('Transmit Time Error (m) ')

    figure
    plot((reference.receiveTime - reference.receiveTime(1,1))',(reference.transmitTime - scalar.transmitTime)')
    title('Transmit Time Error (s) ')

    figure
    plot(nav.cnoL1')
    title('Estimated Carrier To Noise')

    figure
    plot((reference.pos - nav.lsPos)')
    title('Carrier Phase Postion Error')
    legend('X','Y','Z')

    figure
    plot((reference.vel - nav.lsVel)')
    title('Carrier Phase Velocity Error')
    legend('X','Y','Z')

    figure
    plot((reference.clkBias - nav.clkBias)')
    title('Carrier Phase Clock Bias Error')

    figure
    plot((reference.clkDrft - nav.clkDrft)')
    title('Carrier Phase Clock Drift Error')

    figure
    hold on
    plot(nav.rangeRes')
    plot(nav.rangeResVar'.^0.5,'--')
    plot(-nav.rangeResVar'.^0.5,'--')
    title('Kalman Filter Range Residuals')

    figure
    hold on
    plot(nav.rateRes')
    plot(nav.rateResVar'.^0.5,'--')
    plot(-nav.rateResVar'.^0.5,'--')
    title('Kalman Filter Rate Residuals')

end

if plotOldData
    cno = 20;
    bw = 5;

    load(sprintf('scalar_threshold_data/cno%d_codebw2_carrierbw%d_wFeedback_stdx100_1.mat',cno,bw))

    figure
    plot((scalar.receiveTime - scalar.receiveTime(1,1))',scalar.carrierError')
    axis([0 10 -2 2])
    xlabel('Time (s) ')
    ylabel('Phase Error (cyc) ')
    title(sprintf('Bandwidth %d Hz    C/N_0 %d dB-Hz   VelStd 0.01 ',bw,cno))
end