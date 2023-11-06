% we define here a really crude True-FISP a.k.a. bSSFP sequence
% there is no user control for TR/TE, you just specify the ADC time and the
% RF parameters and the rest is calculated to find the fastest posible
% timing. The sequence intensively uses the extended trapezoid
% functionality to achieve near-optimal timing. Due to the requirement for
% splitting the sequence into blocks the TR is increased by approximately
% 40-60 us (rfDeadTime+adcDeadTime) in comparison to the true minimum TR


TR = 13.3e-3;
delta_te = 1.5e-3;
% number of contrasts
n_echos = 5;
TE = zeros(1,n_echos);

% set system limits
% had to slow down ramps and increase adc_duration to avoid stimulation (safe values:13;40;; for plot: 25,110)
sys = mr.opts('MaxGrad',13,'GradUnit','mT/m',...
    'MaxSlew',40,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
    'adcDeadTime', 20e-6);

seq=mr.Sequence(sys);              % Create a new sequence object
fov=300e-3; Nx=32; Ny=32;     % Define FOV and resolution


% ADC duration (controls TR/TE)
% adc_dur=2560; %us
adc_dur=840; %us

% RF parameters 
alpha=60; % deg
thick=20; %mm
rf_dur=1260; % us
rf_apo=0.5;
rf_bwt=2.7;

% Create 'alpha' degree slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',rf_dur*1e-6,...
    'SliceThickness',thick*1e-3,'apodization',rf_apo,'timeBwProduct',rf_bwt,'system',sys);
% consistent with IDEA
gz.amplitude = -gz.amplitude;
gzReph.amplitude = -gzReph.amplitude;

% Define other gradients and Bipolar ADC events
deltak=1/fov;
gxp = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6,'system',sys);
gxn = mr.scaleGrad(gxp,-1);
adc = mr.makeAdc(Nx,'Duration',gxp.flatTime,'Delay',gxp.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',gxp.area/2,'system',sys);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;

% now we have to reshuffle gradients to achieve a half-way optimal timing
% new gz will consist of two parts: 
% 1: slice refocusing from the previous TR followed by slice selection 
%    including the plato an a small bit of the ramp-down
% 2: the remainder of the ramp-down and the slice refocusing for the next TR
gz_parts=mr.splitGradientAt(gz,mr.calcDuration(rf));
gz_parts(1).delay=mr.calcDuration(gzReph);
gz_1=mr.addGradients({gzReph,gz_parts(1)},'system',sys);
[rf,~]=mr.align('right',rf,gz_1);
gz_parts(2).delay=0;
gzReph.delay=mr.calcDuration(gz_parts(2));
gz_2=mr.addGradients({gz_parts(2),gzReph},'system',sys);

% new gr will consist of two parts: 
% 1: prephaser followed by a part of the read gradient including the 
%    beginning of the ramp-down
% 2: the remainer of the ramp-down and the second "prephaser"
gx_parts=mr.splitGradientAt(gxp,ceil(mr.calcDuration(adc)/sys.gradRasterTime)*sys.gradRasterTime);
gx_parts(1).delay=mr.calcDuration(gxPre);
gx_1=mr.addGradients({gxPre,gx_parts(1)},'system',sys);
% adc.delay=adc.delay+mr.calcDuration(gxPre); % we cannot use mr.align here because the adc duration maz be not aligneed to the grad raster
gx_parts(2).delay=0;
gxPre.delay=mr.calcDuration(gx_parts(2));
gx_2=mr.addGradients({gx_parts(2),gxPre},'system',sys);


% Calculate timing
%pe_dur=min(max(tp2(end),tpr1(end-2)),max(tp(end-3),tpr2(end))); 
gxPre.delay=0; % otherwise duration below is misreported
%pe_dur=max(mr.calcDuration(gz_2),mr.calcDuration(gxPre)+gx.riseTime);
% pe_dur=mr.calcDuration(gx_2);
pe_dur = mr.calcDuration(gxPre);

% adjust delays to align objects
% gz_1.delay=max(mr.calcDuration(gx_2)-rf.delay,0);
gz_1.delay=max(mr.calcDuration(gx_2)-rf.delay+rf.ringdownTime,0); % this rf.ringdownTime is needed to center the ADC and the gradient echo in the center of RF-RF period

% rf.delay=rf.delay+gz_1.delay;

% finish timing calculation
% TR=mr.calcDuration(gz_1)+mr.calcDuration(gx_1);
TE((n_echos+1)/2)=TR/2;

for i = 1:n_echos
    if i ~= (n_echos+1)/2
        TE(i) = TE((n_echos+1)/2) + (i - (n_echos+1)/2)*delta_te;
    end
end


% Calculate timing (need to decide on the block structure already)
% helperT is from central of the pulse to end of phase encoding
helperT=ceil((gz.fallTime + gz.flatTime/2 + gzReph.flatTime + gzReph.riseTime ...
        + gzReph.fallTime - mr.calcDuration(gz_2) + max(mr.calcDuration(gz_2),pe_dur) ...
        )/seq.gradRasterTime)*seq.gradRasterTime;

delayTE = zeros(size(TE)) ;

% delayTE(1) is from the end of phase encoding to the beginning of GxPre
% delayTE(>1) is the gap between two Gx
delayTE(1)=TE(1) - helperT - mr.calcDuration(gxPre) - mr.calcDuration(gxp)/2;
for c=2:length(TE)
    delayTE(c) = TE(c) - TE(c-1) - mr.calcDuration(gxp);
end
assert(all(delayTE(1)>=mr.calcDuration(gxPre,gzReph)));
assert(all(delayTE(2:end)>=0));

% delayTR=round((TR - (mr.calcDuration(gz_1) -gz_1.delay) - max(mr.calcDuration(gz_2),pe_dur) ...
%     - (sum(delayTE) - delayTE(1)+2*mr.calcDuration(gxPre)) ...
%     - mr.calcDuration(gxp)*length(TE))/seq.gradRasterTime)*seq.gradRasterTime;

% here delayTR = delayTE(1)
delayTR=round((TR - helperT * 2 - (sum(delayTE) +2*mr.calcDuration(gxPre)) ...
    - mr.calcDuration(gxp)*length(TE))/seq.gradRasterTime)*seq.gradRasterTime;

assert(delayTR>=0);


%% Add Block
% alpha / 2 preparation: shorter TR and half-angle, no PE, no RO
% create 0.5*alpha prep pulse
rf05=rf;
rf05.signal=0.5*rf.signal;

seq.addBlock(rf05,gz_parts(1));
seq.addBlock(gz_2);
% the following delay calculation fails for agressive sequence timing
prepDelay=mr.makeDelay(round((TR/2-mr.calcDuration(gz_1)-mr.calcDuration(gz_2))/ ...
    sys.gradRasterTime)*sys.gradRasterTime); % I know this round() is not 100% correct
seq.addBlock(prepDelay);

% Loop over phase encodes and define sequence blocks
for i=1:Ny
    rf.phaseOffset=pi*mod(i,2);
    adc.phaseOffset=pi*mod(i,2);
        
    gyPre_2 = mr.makeTrapezoid('y','Area',-phaseAreas(i),'Duration',pe_dur,'system',sys); % current PE step
    if i>1
        gyPre_1 = mr.makeTrapezoid('y','Area',phaseAreas(mod(i+Ny-2,Ny)+1),'Duration',pe_dur,'system',sys); % previous PE step
        seq.addBlock(rf,gz_1, gyPre_1);
    else
        seq.addBlock(rf,gz_1);
    end


    % readout prephase
    seq.addBlock(mr.align('left',gyPre_2,gz_2));
    seq.addBlock(mr.makeDelay(delayTE(1)));

    seq.addBlock(gxPre);

    for c=1:length(TE) % loop over TEs
        if mod(c,2)~=0
            gx=gxn;
        else 
            gx=gxp;
        end

        if (c==1)
            seq.addBlock(gx, adc); 
        elseif (c ~= size(TE,2) && c ~= 1)
            seq.addBlock(mr.makeDelay(delayTE(c)));
            seq.addBlock(gx, adc);
        else
            seq.addBlock(mr.makeDelay(delayTE(c)));
            seq.addBlock(gx, adc);
        end

        if c == length(TE)
            seq.addBlock(gxPre)
        end
        % to check/debug TE calculation with seq.testReport() comment out
        % the above line and uncommend the line below this comment block; 
        % change 'c==3' statement to select the echo to test
        %if c==2, seq.addBlock(gx,adc); else, seq.addBlock(gx); end
    end

    seq.addBlock(mr.makeDelay(delayTR));
end
% % finish the x-grad shape 
% seq.addBlock(gxPre);

% check that the calculated TR was reached
% alpha / 2 prep takes 3 blocks
% assert(TR==delayTR + (mr.calcDuration(seq.getBlock(4))+mr.calcDuration(seq.getBlock(5))));

fprintf('Sequence ready\n');
fprintf('TR=%03f ms  TE=%03f ms\n', TR*1e3, TE*1e3);

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% prepare export

seq.setDefinition('FOV', [fov fov thick*1e-3]);
seq.setDefinition('Name', 'mebssfp');

seq.write('mebssfp.seq')       % Write to pulseq file

%seq.install('siemens');

%% plots and displays
seq.plot('timeDisp','us');

%% plot entire interpolated wave forms -- good for debugging of gaps, jumps,
% etc, but is relatively slow
%gw=seq.gradient_waveforms();
%figure; plot(gw'); % plot the entire gradient waveform
gw_data=seq.waveforms_and_times();
figure; plot(gw_data{1}(1,:),gw_data{1}(2,:),gw_data{2}(1,:),gw_data{2}(2,:),gw_data{3}(1,:),gw_data{3}(2,:)); % plot the entire gradient waveform
title('gradient waveforms');

%% trajectory calculation 
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces

figure; plot(t_ktraj,ktraj'); title('k-space components as functions of time'); % plot the entire k-space trajectory
figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
title('2D k-space');

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  

rep = seq.testReport;
fprintf([rep{:}]);

