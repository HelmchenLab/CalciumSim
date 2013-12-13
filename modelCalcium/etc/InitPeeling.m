function [ca_p,exp_p,peel_p, data] = InitPeeling(dff, rate)

% ca_p: parameters of elementary (1 AP) calcium transient
ca_p.onsetposition =0.0;    % onset position(s)
ca_p.onsettau=0.0000001;    % onset tau (s)
ca_p.offset=0;              % baseline offset (%)
ca_p.amp1=2.5;              % amplitude 1  (%)
ca_p.tau1=0.6;              % tau1 (s)
ca_p.amp2=0;                % amplitude 2 (%)
ca_p.tau2=1.0;              % tau2 (s)
ca_p.integral=0.0;          % integral below curve (%s)
ca_p.scale=1.0;             % scale factor to scale entire trace (s)

% exp_p: experiment parameters, including dye properties and data acquisition 
exp_p.numpnts = length(dff); % numpoints
exp_p.acqrate = rate;        % acquisition rate (Hz)
exp_p.noiseSD = 1.2;        % noise stdev of DF/F trace (in percent), should be specified by the user
exp_p.indicator = 'OGB-1';  % calcium indicator
exp_p.dffmax = 93;          % saturating dff max (in percent)
exp_p.kd = 200;             % dye dissociation constant (nM)
exp_p.carest = 50;          % presumed resting calcium concentration (nM)

% peel_p: parameters for peeling algorithm
peel_p.padding = 20;        % number of points for padding before and after
% peel_p.sdnoise = 1.4;       % expected SD baseline noise level
peel_p.smtthigh = 2.4;      % Schmitt trigger - high threshold (multiple of exp_p.noiseSD)
peel_p.smttlow = -1.2;      % Schmitt trigger - low threshold (multiple of exp_p.noiseSD)
peel_p.smttbox= 3;          % Schmitt trigger - smoothing box size (in points)
peel_p.smttmindur= 0.3;     % Schmitt trigger - minimum duration (s)

% HL: 2012-05-04
% new parameter: max. frames fro smttmindur
% if the frame rate is high, number of frames for smttmindur can be
% large, thereby increasing false negatives
% if smttminFrames is set, use binning to reduce the number of
% frames to this value for high frame rates
% peel_p.smttminFrames = 20;

peel_p.smttnumevts= 0;      % Schmitt trigger - number of found events
peel_p.slidwinsiz= 10.0;    % sliding window size - event detection (s)
peel_p.maxbaseslope= 0.5;   % maximum baseslope %/s
peel_p.evtfound=0;          % flag - 1: crossing found 
peel_p.nextevt=0;           % next crossing found (s)
peel_p.nextevtframe=0;      % next crossing found (frame number)
peel_p.intcheckwin=0.5;     % window to the right - for integral comparison (s)
peel_p.intacclevel=0.5;     % event integral acceptance level (0.5 means 50%)
peel_p.fitonset=0;          % flag - 1: do onset fit, only useful if 1/frameRate <= rise of CacliumTransient
peel_p.fitwinleft=0.5;     % left window for onset fit (s)
peel_p.fitwinright=0.5;    % right window for onset fit (s)
peel_p.negintwin=0.5;       % window to the right - for negativeintegral check(s)
peel_p.negintacc=0.5;       % negative acceptance level (0.5 means 50%)
peel_p.stepback=1.5;        % stepsize backwards for next iteration (s)
peel_p.fitupdatetime=2;     % how often the linear fit is updated (s)

% data: data struct 
data.dff = dff;
data.tim = 1:length(data.dff); 
data.tim = data.tim./exp_p.acqrate;
data.intdff = 1:length(data.dff);                % integral curve
data.singleTransient = zeros(1,exp_p.numpnts);
data.model = zeros(1,exp_p.numpnts);
data.spiketrain = zeros(1,exp_p.numpnts);
data.slide = zeros(1,exp_p.numpnts);            % sliding curve, zero corrected
data.temp = 1:length(data.dff);                 % temporary wave
data.peel = zeros(1,exp_p.numpnts);
data.peel = data.dff;
data.spikes = zeros(1,1000);                    % array for found spikes times
data.numspikes = 0;                             % number of spikes found


