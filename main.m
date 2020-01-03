
clear all, close all

basepath = '/home/mario/MEGA/Neuroscience/';
addpath([basepath 'Awareness_Perception_Recognition/Metacognition/pjE3.47-AccessingAwareness/'])
addpath([basepath 'MotorControl_Agency_FreeWill/Pupillometry_EyeTracking/'])
addpath([basepath 'Awareness_Perception_Recognition/Psychophysics_SignalDetectionTheory/'])
datapath = [basepath 'Data/Pupillometry/'];	

toplt_smoodat = 0; 
toplt_evlkapd = 0;
toplt_pcbiplot = 0;
todobatch47 = 4; % 0 4 7

global fs decf pxPmm 
pxPmm = 10; % px to mm ratio
fs = 240; % sampling rate in Hz
nyl = fs / 2; % Nyquist limit 
decf = 1; % downsampling rate
fs2 = fs / decf;

% window parameters
lkev = {'onset_fp', 'onset_st', 'offset_st', 'onset_qu1', 'rt1', 'onset_qu2', 'rt2'};
blp = 100; % baseline period in ms 
bls = blp * fs2 / 1000; % baseline period in samples: 24 for 100ms
% 0 to 2000
rwp2k = 2000; % response window length in ms (prev. 3000ms) 
rws2k = rwp2k * fs2 / 1000; % response window length in samples: 600 for 2500ms
rwt2k = (1:rws2k) * 1000 / fs2;
% 500 and 1200 based on Hoeks & Levelt pupillary response function
awp1 = 500; awp2 = 1100; % analysis window start and end in ms  
aws1 = awp1 * fs2 / 1000; aws2 = awp2 * fs2 / 1000; % analysis window start and end in sample number
awt = (aws1:aws2) * 1000 / fs2;
% time window where pirf is more than 20% of maximum
a20s = find(pirf(rwt2k) > max(pirf(rwt2k)) * 0.2); 
a20s1 = a20s(1); a20s2 = a20s(end);
a20t = a20s * 1000 / fs2;
% chosen analysis window
aws = a20s;
awt = a20t;

switch todobatch47
case 0, ss = []; ms = [1 3]; is = 1:3; %ss = randi(10); ms = 1+2*(randi(2)-1); is = randi(3);
case 4, ss = 1:6; ms = [1 3]; is = 1:3; dsf = {'ftr1', 'ftr2', 'qu1', 'qu2', 'qu12'};
case 7, ss = 7:10; ms = [1 3]; is = 1:3; dsf = {'ftr', 'qu1', 'qu2', 'qu12', 'noi'};
end

load tables47.mat % Ft

%% read data
readRawData
% PA{s,m}, MA{s,m}, i_eventsA{s,m}
% DS{s,m,i}, T{s,m,i}, Fm{s,m,i}, Ft_smi{s,m,i}


%% explore raw data
plt_raw(PA{s, m}, MA{s, m}, i_eventsA{s, m}, s, i);


%% plot response window waves and pool data
plotAndorPoolPupilAmp


%% pool data into tables
poolIntoTables         
% Ft_mi = {}; Fma = {}; rti_qu1a = {}; rti_qu2a = {};  


%% completing Ft table 
completeFtTable
% Ft


%% plot power spectral density
plt_spectra(T{s,m,i}, DS{s,m,i});


%% event window length averages: fpdi;(2) sdi; onset_qu1 - offset_st;(4) rti1; onset_qu2 - rt1;(6) rti2
[iei, avg_iei, std_iei] = plt_evwist(timelog{s,m,i}, 1);


%% comparing waveforms 1 
%  printing sample-by-sample massive t-testing
%  fit loads of regression models and anovas 
%  the proper way would be cluster-based permutation test
testingWaveformDifferences


%% comparing waveforms 2
%  plotting within 500-1100 ms window post event, individual averages 
%  analyze lkev{2,4,5,6,7}, and dsf{3,5} for each subject
%  calculate mean amplitude in awt: absolute values and differences
close all

ss = 1:6; m = 1; i = 1; ei = 2; fi = 3; lkev{ei}, dsf{fi}

if s < 7, [facvec, levs] = make_facvec(Fm{s,m,i}(:,fi), s);
else,     [facvec, levs] = make_facvec(Fm{s,m,i}(:,fi), s, Fm1{s,m,i}(:,fi));
end

% plot subject average pupil size waves
for s = ss
  if s < 7
    plt_evlk_apd(Dav_smief{s,m,i,ei,fi}, Dse_smief{s,m,i,ei,fi}, rwt2k, ...
        facvec, lkev{ei}, dsf{fi}, 1)
    A = cellfun(@(x) mean(x(a20s)), Dav_smief{s,m,i,ei,fi}); diff(A)
  else
    plt_evlk_apd(Dav_smief{s,m,i,ei,fi}, Dse_smief{s,m,i,ei,fi}, rwt2k, ...
        [facvecc{1}; facvecc{2}], lkev{ei}, dsf{fi}, 1)
  end
end


%% Metacognition analysis

close all
Ft = Ft4(Ft4.ftr2~=categorical(2), :);
is4 = 1;
seom = @(x) std(x, 'omitnan') ./ sqrt(size(x, 1));
cr2vr = @(x) [sum(x(1:4)) 2/3*sum(x(5:6)) 2/3*sum(x(6:7)) 2/3*sum(x(7:8))];

makeTableG % Subject-averaged table G

% CHOOSE
subgrp = 'avg' % one, avg, pool, all
task = 3       % 1:3
mj = 3        % [1 3]  
islr = 1       % 0 1 % S1,S2 correspond to L,R? 
rstr = '_rS1'      % response specific?: '' '_rS1' '_rS2'  
ist1uv = 0   % type1 unequal variance model?
ist2uv = 1, mj_s = 3  % type2 unequal variance model? NEEDS BEFORE TO RUN t1uv. 1: CR-based, 3: VR-based
% don't modify these (usually):
s1divs2=1; iscr2vr=0; iseqvar=0; nR_thr=4; AUC={}; ncols=10; nrows=3;

switch subgrp
case 'one',  si = ss(1);
case 'avg',  si = ss; sa = ss(end)+1;   
case 'pool', si = ss;
case 'all',  si = ss;
end


for tk = task
    
  switch  tk
    case 1, lvs = 2:8; if mj==3 || iscr2vr, nR_thr = 1; end
    case 2, lvs = 2:4; %dpc lvs = 6:8; S1 = [4 5]; S2 = @(x) [x 9-x];
    case 3, lvs = 2:3; 
  end
  
  % Fit 
  clear fit
  
  % fetch perceptual decisions
  % IMPORTANT: qu1 is YN for detec, but accuracy for discr1 and discr2
  q1 = Ft{any(Ft.sid==si, 2) & Ft.task==tk & Ft.mjdg==mj, 'qu1'};  
  q2 = Ft{any(Ft.sid==si, 2) & Ft.task==tk & Ft.mjdg==mj, 'qu2'};  
  q12 = Ft{any(Ft.sid==si, 2) & Ft.task==tk &   Ft.mjdg==mj, 'qu12'};      
  %uq1 = unique(rmmissing(q1))'; uq2 = unique(rmmissing(q2))'; uq12 = unique(rmmissing(q12))';
  uq1 = [0 1]; uq2 = 1:4; uq12 = [-4:-1 1:4];
  
  %Fit meta-d' model
  fitMetadModel  % fit
  
  
  % meta-d' plots
  % fig04a, fig05a, fig08ac, figC3ac(iseqvar=1), figC1a(iseqvar=1), figC4abcd, figC5abcde
  subplot(nrows, ncols, 9 + ncols*(tk-1))   
  switch tk
      case 1, axis([lvs(1)-0.2 lvs(end)+0.2 -2 3.5])
        if ~isempty(rstr), ylim([-3.5 5]), end                        
      case 2, axis([lvs(1)-0.2 lvs(end)+0.2 -2.5 4])
        if ~isempty(rstr), ylim([-3.5 5]), end 
  end, hold on, grid on, box on
  if (mj ~= 3 || tk ~= 1) && ~ist1uv
    lgs = {'d_a', ' meta-d_a', ' meta-d_a - d_a'}; pl = [];
    if strcmp(subgrp, 'avg')
      eb1 = errorbar(lvs, [fit{sa}.da], [fit{sa+1}.da], 'k');
      eb2 = errorbar(lvs, getfitfields(fit{sa}, ['meta_da' rstr], lvs), getfitfields(fit{sa+1}, ['meta_da' rstr], lvs), '-'); 
      eb3 = errorbar(lvs, getfitfields(fit{sa}, ['M_diff' rstr], lvs), getfitfields(fit{sa+1}, ['M_diff' rstr], lvs), 'k--');   
      pl = [eb1 eb2 eb3];
    elseif strcmp(subgrp, 'one') || strcmp(subgrp, 'pool')
      lvs_p = lvs;
      if length([fit{1}.da]) ~= length(lvs), lvs_p(1) = []; end
      pl = plot(lvs_p, [fit{1}.da], lvs_p, getfitfields(fit{1}, ['meta_da' rstr], lvs_p), '--', lvs_p, getfitfields(fit{1}, ['M_diff' rstr], lvs_p), ':');       
    elseif strcmp(subgrp, 'all') 
      for s = ss 
        lvs_p = lvs;
        if length([fit{s}.da]) ~= length(lvs), lvs_p(1) = []; end
        pl(3*s-2:3*s) = plot(lvs_p, [fit{s}.da], lvs_p, getfitfields(fit{s}, ['meta_da' rstr], lvs_p), '--', lvs_p, getfitfields(fit{s}, ['M_diff' rstr], lvs_p), ':');
        hold on
        lgs{3*s-2} = ['s' num2str(s) ' d_a'];       
        lgs{3*s-1} = ['s' num2str(s) ' meta-d_a'];       
        lgs{3*s} = ['s' num2str(s) ' meta-d_a - d_a'];       
      end      
    end
    %if ~strcmp(subgrp,' all') && tk==3, legend(pl, lgs, 'Location', 'best'), end
    title({'Discriminability d_a and meta-d_a'})
    lgs = {'d_a' 'meta-d_a' 'meta-d_a - d_a'}; 
    ca = gca; 
    set(ca.Children(2), 'Color', [.5 .5 .5])
  elseif ist1uv % figB3a, figB2a, figB4a, figB4c
    lgs = {'d_a', 'd_1', 's'}; pl = []; 
    switch tk
      case 1, ylim([-2 5])
      case 2, ylim([-2.5 4.5])
    end
    if strcmp(subgrp, 'avg')      
      eb1 = errorbar(lvs, [fit{sa}.da], [fit{sa+1}.da], 'k');
      eb2 = errorbar(lvs, [fit{sa}.d_1], [fit{sa+1}.d_1], 'k-.');
      eb3 = errorbar(lvs, [fit{sa}.s], [fit{sa+1}.s], 'k:'); 
      pl = [eb1 eb2 eb3];
    end
    title({'Discriminabilities d_a and d_1'', and s.d.r.'})
    ca = gca;     
  end
  for c = ca.Children, set(c, 'LineWidth', 2), end
  xlabel('Signal strength'), ylabel('s.d. units')
  legend(pl, lgs, 'Location', 'best')
  
  
  % fig04b, fig05b, fig08bd, figC3bd(iseqvar=1), figC1b(iseqvar=1)
  subplot(nrows, ncols, 10 + ncols*(tk-1))
  switch tk
      case 1, axis([lvs(1)-0.2 lvs(end)+0.2 -2 3.5])
      case 2, axis([lvs(1)-0.2 lvs(end)+0.2 -2.5 4])
  end, hold on, grid on, box on
  if (mj ~= 3 || tk ~= 1) && ~ist1uv
    if strcmp(subgrp, 'avg') 
      eb1 = errorbar(lvs, getfitfields(fit{sa}, [t1ORmeta rstr], lvs), getfitfields(fit{sa+1}, [t1ORmeta rstr], lvs), 'k-'); 
      eb2 = errorbar(repmat(lvs', [1 3]), reshape([fit{sa}.t2ca_rS1], 3, length(lvs))', reshape([fit{sa+1}.t2ca_rS1], 3, length(lvs))', 'k--'); 
      eb3 = errorbar(repmat(lvs', [1 3]), reshape([fit{sa}.t2ca_rS2], 3, length(lvs))', reshape([fit{sa+1}.t2ca_rS2], 3, length(lvs))', 'k--'); 
      pl = [eb1 eb2 eb3];
    else
      for s = ss
        lvs_p = lvs;
        if length(getfitfields(fit{s}, [t1ORmeta rstr], lvs)) ~= length(lvs), lvs_p(1) = []; end  
        pl = plot(lvs_p, getfitfields(fit{s}, [t1ORmeta rstr], lvs), '-*', ...
          lvs_p, reshape([fit{s}.t2ca_rS1], 3, length(lvs_p))', '-.', ...
          lvs_p, reshape([fit{s}.t2ca_rS2], 3, length(lvs_p))', '--'); hold on
        if strcmp(subgrp, 'one') || strcmp(subgrp, 'pool'), break, end 
      end
    end
    %if tk==3, legend(pl, lgs, 'Location', 'best'), grid on, end
    title({'Type 1 and 2 criteria for meta-d_a fit'})
    lgs = {'type2c_a' 'meta-c_a'};
  elseif ist1uv  %figB3b, figB2b figB4b, figB4d
    lgs = {'c_1'}; 
    switch tk
      case 1, ylim([-2 5])
      case 2, ylim([-2.5 4.5])
    end
    if strcmp(subgrp, 'avg')      
      errorbar(lvs, [fit{sa}.c_1], [fit{sa+1}.c_1], 'k-*');     
      nthres = length(fit{sa}(end).c_1_all);
      pl = errorbar(repmat(lvs', [1 nthres]), reshape([fit{sa}.c_1_all], nthres, length(lvs))', reshape([fit{sa+1}.c_1_all], nthres, length(lvs))', 'k--');       
    end
    title({'Type 1 criteria'}), grid on
  end
  ca = gca; 
  for c = ca.Children, set(c, 'LineWidth', 2), end
  xlabel('Signal strength'), ylabel('s.d. units')
  %legend(pl, lgs, 'Location', 'northeast')
  
   
  % stimulus strength levels
  calcSignalwiseVars;
  %q1={}; q2={}; q12={}; idx1={}; jdx1={}; idx2={}; jdx2={}; idx12={}; jdx12={};
  %labels1={}; labels2={}; labels12={}; labels32={}; scores1={}; scores2={}; scores12={};  
  %nwq1 = linspace(0, 1, 2); nwq12 = linspace(0, 1, 8); nwq2 = linspace(0, 1, 4);
   
  % fig02ad: Psychometric curve: hits and falsal
  if tk==1, subplot(nrows, ncols, 1 + ncols*(tk-1))
    ylabel('Accuracy'), xlabel('Signal strength')
    axis([0.8 lvs(end)+0.2 0 1]), grid on, box on, hold on
    if mj == 1, title('Psychometric 2-points visibility')
      switch subgrp
        case 'all',  plot([1 lvs], q1M(ss,[1 lvs]))   
        case 'avg',  errorbar([1 lvs], mean(q1M), seom(q1M))
        case 'pool', errorbar([1 lvs], q1mp, q1sp)
        case 'one',  plot([1 lvs], q1M(si,[1 lvs]))   
      end
    else,       title('Psychometric 2-points visibility')
      switch subgrp
        case 'all',  plot([1 lvs], q21M(ss,[1 lvs]))   
        case 'avg',  errorbar([1 lvs], mean(q21M), seom(q21M))
        case 'pool', errorbar([1 lvs], q21mp, q21sp)
        case 'one',  plot([1 lvs], q21M(si,[1 lvs]))   
      end    
    end
    line([1 lvs], .5*ones(1, length(lvs)+1))
    ca = gca; 
    set(ca.Children(1), 'Color', [.5 .5 .5], 'LineStyle', ':', 'LineWidth', 1) 
    set(ca.Children(2), 'LineWidth', 2, 'Color', 'k')
    %print('deliveries/fig1ad.eps', '-deps')
  end
  

  % fig02be, fig06abde, fig09ac: Psychometric curve: 2nd order judgment
  subplot(nrows, ncols, 2 + ncols*(tk-1))
  axis([0.8 lvs(end)+0.2 0 1]), grid on, box on, hold on
  switch mj 
    case 1, title('Psychometric 4-points confidence'), ylabel('Confidence (CR)')      
    case 3, title('Psychometric 4-points visibility'), ylabel('Visibility (FVR)'),
  end, xlabel('Signal strength')
  switch subgrp
    case 'all', plot([1 lvs], q2M(ss,[1 lvs])), hold on   
      if mj==1||tk~=1, plot([1 lvs], q2M0(ss,[1 lvs])), plot([1 lvs], q2M1(ss,[1 lvs])), end
    case 'avg', errorbar([1 lvs], mean(q2M, 'omitnan'), seom(q2M)), hold on
      if mj==1||tk~=1, errorbar([1 lvs], mean(q2M0, 'omitnan'), seom(q2M0)), errorbar([1 lvs], mean(q2M1, 'omitnan'), seom(q2M1)), end    
    case 'pool',errorbar([1 lvs], q2mp, q2sp), hold on
      if mj==1||tk~=1, errorbar([1 lvs], q2mp0, q2sp0), errorbar([1 lvs], q2mp1, q2sp1), end
    case 'one', plot([1 lvs], q2M(si,[1 lvs])), hold on   
      if mj==1||tk~=1, plot([1 lvs], q2M0(si,[1 lvs])), plot([1 lvs], q2M1(si,[1 lvs])), end 
  end 
  line([1 lvs], .5*ones(1, length(lvs)+1))
  ca = gca; 
  switch tk
    case 1
      if mj==1, legend(ca, 'Yes-CR', 'No-CR', 'CR', 'Location', 'southeast'), end
      set(ca.Children(1), 'Color', [.5 .5 .5], 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility', 'off') 
      set(ca.Children(1), 'LineWidth', 2, 'Color', 'k')
      set(ca.Children(2), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k')
      set(ca.Children(3), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k')
    case 2
      set(ca.Children(1), 'Color', [.5 .5 .5], 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility', 'off')             
      set(ca.Children(1), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k')
      set(ca.Children(2), 'LineWidth', 2, 'LineStyle', '-.', 'Color', 'k')
      set(ca.Children(3), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'k') 
      switch islr   
        case 1, legend(ca, 'All', 'CCW', 'CW', 'Location', 'southeast')
        case 0, legend(ca, 'All', 'Incorrect', 'Correct', 'Location', 'southeast')  
            set(ca.Children(2), 'LineStyle', ':')
      end  
    case 3    
      legend(ca, 'All', 'Incorrect', 'Correct', 'Location', 'northwest')  
      set(ca.Children(1), 'Color', [.5 .5 .5], 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility', 'off')       
      set(ca.Children(1), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k')
      set(ca.Children(2), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k')
      set(ca.Children(3), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'k')     
  end

  
  % figA1abcdef
  axis([0.9 lvs(end)+0.1 0 1]), grid on, box on, hold on
  for i=1:4, eb(i)=errorbar([1 lvs], mean(q2Mrating{i}, 'omitnan'), seom(q2Mrating{i})); end % q2M0rating  q2M1rating
  legend({'1', '2', '3', '4'}, 'Location', 'best')
  ca = gca;
  for c=ca.Children, set(c, 'LineWidth', 2), end
  %scatter(repmat(1:8,1,6),reshape(log10(q2M0)',1,[])), lsline
  %scatter(repmat(1:8,1,6),reshape(log10(q2M1)',1,[])), lsline
  %fl=fitlm(repmat(1:8,1,6),reshape(log10(q2M0)',1,[])); fl.Rsquared

 
  % fig02cf, fig06cf, fig09bd: Psychometric accuracy (type I) & Metapsychometric accuracy (type II)
  subplot(nrows, ncols, 3 + ncols*(tk-1))   
  if tk~=3, axis([lvs(1)-0.2 lvs(end)+0.2 0.4 1]), grid on, box on, hold on
  else,     axis([lvs(1)-0.2 lvs(end)+0.2 0.1 1]), grid on, box on, hold on
            plot(lvs, .25*ones(1, length(lvs)), ':')
  end
  CF = cellfun(@mean, labels2);
  switch subgrp
  case 'all', plot(lvs, CF(ss,lvs)) 
  case 'avg', errorbar(lvs, mean(CF(ss,lvs)), seom(CF(ss,lvs)))
  case 'pool', CFsl = cellfun(seom, labels2); errorbar(lvs, CF(lvs), CFsl(lvs))
  case 'one', plot(lvs, CF(si,lvs))   
  end
  ylabel('Accuracy'), xlabel('Signal strength')
  if mj==1, title('Task and introspective accuracy')
  else,     title('Task accuracy'), end
  line(lvs, .5*ones(1, length(lvs)))
  ca = gca; 
  if tk ~= 3
    set(ca.Children(1), 'Color', [.5 .5 .5], 'LineStyle', ':', 'LineWidth', 1.5, 'HandleVisibility', 'off') 
    set(ca.Children(1), 'LineWidth', 2, 'Color', 'k')
  else
    set(ca.Children(1), 'Color', [.5 .5 .5], 'LineStyle', ':', 'LineWidth', 1.5, 'HandleVisibility', 'off') 
    set(ca.Children(2), 'Color', [.5 .5 .5], 'LineStyle', ':', 'LineWidth', 1.5, 'HandleVisibility', 'off') 
    set(ca.Children(1), 'LineWidth', 2, 'Color', 'k')
  end
  if mj==1 || tk~=1
    CF = cellfun(@mean, labels32);
    switch subgrp
    case 'all', plot(lvs, CF(ss,lvs))   
    case 'avg', errorbar(lvs, mean(CF(ss,lvs)), seom(CF(ss,lvs)))       
    case 'pool', CFsl = cellfun(seom, labels32); errorbar(lvs, CF(lvs), CFsl(lvs)) 
    case 'one', plot(lvs, CF(si,lvs))   
    end    
    set(ca.Children(1), 'Color', [.5 .5 .5], 'LineWidth', 2) 
    legend(ca, 'Task accuracy', 'Introspective accuracy', 'Location', 'northwest')
  end
  %if mj==1 || tk~=1, subplot(nrows, ncols, 4 + ncols*(tk-1))   
  %  axis([lvs(1)-0.2 lvs(end)+0.2 0.3 1.1]), grid on, hold on
  %  plot(lvs, .5*ones(1, length(lvs)), ':')
  %  CF = cellfun(@mean, labels32);
  %  switch subgrp
  %  case 'all', plot(lvs, CF(ss,lvs))   
  %  case 'avg', errorbar(lvs, mean(CF(ss,lvs)), seom(CF(ss,lvs)))       
  %  case 'pool', CFsl = cellfun(seom, labels32); errorbar(lvs, CF(lvs), CFsl(lvs)) 
  %  case 'one', plot(lvs, CF(si,lvs))   
  %  end    
  %end
    
  
  % ROC curves datapoints lie in the same abscissa because we reuse the
  % same reference for null stimulus (level 1)
  for rty = 1:4 % ROC type
    if (mj==3 && tk==1 && (rty==3 || rty==4)) || (tk~=1 && rty==4), continue, end
    subplot(nrows, ncols, (rty+4)+ncols*(tk-1))
            
    lgs = {}; eb = []; 
    for l = lvs
      X={}; Y={}; T={}; 
      for s = ss
        switch subgrp
        case {'pool', 'one'}, NBootVal = 1000;
        otherwise,            NBootVal = 0; 
        end
        
        switch rty
        case 1  % Type 1 sensitivity
          % here, the score threshold is referenced to the FPR (chose yes when no)    
          [X{s} Y{s} T{s} AUC{rty,tk,l,s} OPTOP NCY NCYN] = perfcurve(...
              labels1{s,l}, scores1{s,l}, 1, 'NBoot', NBootVal, 'TVals', 'all');
          if strcmp(subgrp, 'pool') || strcmp(subgrp, 'one'), auc = AUC{rty,tk,l,s}(2);
          else,                                               auc = AUC{rty,tk,l,s}; end
          G{G.sid==s & G.task==tk & G.mjdg==mj & G.ftr1==l,'auc1'} = auc;
          nsc = 2;
          misthr = setdiff(1:nsc, unique(idx1{s,l}));          
        case 2  % figB1abcd: Pseudo-type 1 sensitivity (type 1 augmented with metadecisions)
          % here, the score threshold is referenced to a 8-point score combining 1 and 2 decisions when no
          [X{s} Y{s} T{s} AUC{rty,tk,l,s} OPTOP NCY NCYN] = perfcurve(...
              labels1{s,l}, scores12{s,l}, 1, 'NBoot', NBootVal, 'TVals', 'all');
          if strcmp(subgrp, 'pool') || strcmp(subgrp, 'one'), auc = AUC{rty,tk,l,s}(2);
          else,                                               auc = AUC{rty,tk,l,s}; end
          G{G.sid==s & G.task==tk & G.mjdg==mj & G.ftr1==l,'auc12'} = auc;
          nsc = 8;
          misthr = setdiff(1:nsc, unique(idx12{s,l}));
        case {3 4}  % figC2abc: Absolute type 2 sensitivity
          % here, the score threshold is referenced to the 2nd order FPR (CR when incorrect) 
          if rty==4 % compute SDI AUC 
            labels2{s,l} = labels2{s,l}(~q1{s,l}); 
            scores2{s,l} = scores2{s,l}(~q1{s,l});
            
            if length(unique(labels2{s,l})) == 0
              continue
              %warning('Skipping: no negative predictions in labels2 s=%u, l=%u', s, l)
              %labels2{s,l} = [0 1 0 1]';              
            elseif length(unique(labels2{s,l})) == 1
              continue
              %warning('Flipping one trial in s=%u, l=%u labels2 for plotting', s, l)
              %labels2{s,l}(randi(length(labels2{s,l}))) = 1 - labels2{s,l}(randi(length(labels2{s,l})));
            end
            if length(unique(scores2{s,l})) == 0
              continue
              %warning('Skipping: no negative predictions in scores2 s=%u, l=%u', s, l)
              %scores2{s,l} = nwq2';
            elseif length(unique(scores2{s,l})) < 4
              continue
              %warning('Flipping one trial in s=%u, l=%u scores2 for plotting', s, l)
              %scores2{s,l}(randi(length(labels2{s,l}),1,4)) = nwq2';
            end
          end        
          
          if all(scores2{s,l}~=0), scores2{s,l}(randi(length(scores2{s,l})))=0; end % patch 
          try
            [X{s} Y{s} T{s} AUC{rty,tk,l,s} OPTOP NCY NCYN] = perfcurve(...
                labels2{s,l}, scores2{s,l}, 1, 'NBoot', NBootVal, 'TVals', 'all');
            if strcmp(subgrp, 'pool') || strcmp(subgrp, 'one'), auc = AUC{rty,tk,l,s}(2);
            else,                                               auc = AUC{rty,tk,l,s}; end
            if rty==3, G{G.sid==s & G.task==tk & G.mjdg==mj & G.ftr1==l,'auc2'} = auc;
            else       G{G.sid==s & G.task==tk & G.mjdg==mj & G.ftr1==l,'auc2sid'} = auc; 
            end
            nsc = 4;
            misthr = setdiff(1:nsc, unique(idx2{s,l}));
          catch ME
            warning('Ignoring exception'), warning(ME.identifier, ME.message)        
          end          
        end
        
        if misthr
          [~, ind] = sort([flip(setdiff(1:nsc+1, misthr)) misthr], 'descend');
          tmp = [X{s}' NaN(size(X{s},2), length(misthr))]'; X{s} = tmp(ind,:);
          tmp = [Y{s}' NaN(size(Y{s},2), length(misthr))]'; Y{s} = tmp(ind,:);
          tmp = [T{s}' NaN(size(T{s},2), length(misthr))]'; T{s} = tmp(ind,:);          
        end
               
        switch subgrp
        case {'pool' 'one'}
          eb(l-1) = errorbar(X{s}(:,1), Y{s}(:,1), Y{s}(:,1)-Y{s}(:,2), Y{s}(:,3)-Y{s}(:,1), '-o'); 
          hold on, break
        case 'all'
          eb(l-1) = plot(X{s}(:,1), Y{s}(:,1), '-o'); 
          pl = plot(0:.01:1, 0:.01:1, ':', X{s}(:,1), T{s}(:,1), '--');     
          hold on
        end
      end       
      
      try
        if strcmp(subgrp, 'avg') 
          [X{sa} X{sa+1}] = grpstats([X{:}]', [], {'mean', 'sem'});
          [Y{sa} Y{sa+1}] = grpstats([Y{:}]', [], {'mean', 'sem'});
          T{sa} = grpstats([T{:}]', [], 'mean');
          [AUC{rty,tk,l,sa} AUC{rty,tk,l,sa+1}] = grpstats([AUC{:}], [], {'mean', 'sem'});
          eb(l-1) = errorbar(X{sa}, Y{sa}, Y{sa+1}); hold  on
          set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(eb(l-1))))
          eb(l-1) = errorbar(X{sa}, Y{sa}, X{sa+1}, 'horizontal'); 
          %pl = plot(0:.01:1, 0:.01:1, ':', X{sa}, T{sa}, ':');               
        else
          pl = plot(0:.01:1, 0:.01:1, ':', X{1}(:,1), T{1}(:,1), ':');       
        end      
        lgs{l-1} = num2str(l);      
      catch ME
        warning('Ignoring exception'), warning(ME.identifier, ME.message)        
      end      
    end 
     
    xlabel('1 - Specificity'), ylabel('Sensitivity'), grid on
    tls = 'thr'; % 'Threshold on classifier scores';
    for i = 1:length(lgs), lgs{i} = ['L' lgs{i}]; end
    legend(eb, lgs, 'Location', 'southeast')  %legend([eb pl(2)], [lgs {tls}], 'Location', 'southeast')   
    switch rty
      case 1, title('Type 1 ROC')
      case 2, title('Pseudo-type 1 ROC') %print('deliveries/figB1.eps', '-deps')
      case 3, title('Absolute Type 2 ROC') %print('deliveries/figC2.eps', '-deps')
      case 4, title({'Absolute','type 2 ROC (SDI)'})
    end    
    ca = gca; cf = gcf; 
    if mj==1, for c = ca.Children, set(c, 'LineWidth', 1.5), end
    else,     for c = cf.Children(2).Children, set(c, 'LineWidth', 1.5), end, set(cf.Children(2), 'YTick', 0:.2:1)
    end
    
  end           
end


%% Confusion matrices

disp({'Hits' 'False alarms'; 'Misses' 'Correct rejections'})
disp({'Rightly sure' 'Overconfident'; 'Underconfident' 'Rightly hesitant'})
for i = is
  % objective sensitivity (Type I)
  T = Ft(Ft.mjdg==1 & Ft.task==i,:);
  ct = crosstab(T.qu1, T.ftr1 > 1); % decision, condition
  cm1 = fliplr(flipud(ct))
  confmvars(cm1)
  % subjective sensitivity (Type II)
  % SDI is negpredval for Type II measures
  T = Ft(Ft.mjdg==1 & Ft.task==i,:);
  ct = crosstab(T.qu2 > 2, (T.ftr1>1 & T.qu1) | (T.ftr1==1 & ~T.qu1)); % metadecision, correctness
  cm2 = fliplr(flipud(ct))
  confmvars(cm2)  
end


%% adding to Ft

Ft = Ft4(Ft4.ftr2~=categorical(2),:);
%Ft.mjdg = double(Ft.mjdg); Ft.task = double(Ft.task);
Ft{Ft.rti1<0,'rti1'}=NaN; % Ft{Ft.rti1<0|Ft.rti1>4,'rti1'}=NaN;
Ft{Ft.rti2<0.0071,'rti2'}=NaN;
Ft.right1 = double(Ft.right1);
Ft.right2 = double(Ft.right2);
Ft.tk2qu1lr = ((Ft.ftr1_raw>4 & Ft.qu1) | (Ft.ftr1_raw<5 & ~Ft.qu1));
%writetable('Ft', 'Ft4f13.csv')

% removing outliers (beyond 3-sigma)
rt4sigma = std(Ft{~any(isnan(Ft{:, [4 5]}), 2), [4 5]}, 1) * 4; 
rt4sigma_avg = mean(rt4sigma);
rt4s_out = any(Ft{:,[4 5]} > rt4sigma_avg, 2); % 268 / 10445 = 2.57%


%% regression 
fittingGlmes


%% a bit of multivariate analysis
[U,S,V] = svd(Ft{Ft.task==1 & Ft.mjdg==1, [3:6 8:10 13 15:18]});


%% Anova
doingAnovas


%% sundry plots

lkev, dsf{fi}

% PS distributions
for ei = 1:eis
  subplot(eis(end),1,ei)
  hist(Ft{:,cnv+ei},100); 
end

%
figure
parallelcoords(Ft(:, [2:5 7:8]), 'group', Ft.qu1)
andrewsplot(Ft(:, [2:5 7:8]), 'group', Ft.qu1)

% fig03, fig07, fig10: errorbar plots of RTs
xmin = [0.5 0.8 0.5]; ymin = [.1 .1 .1];
xmax = [8.5 4.2 3.5]; ymax = [1.4 0.9 1.4];
for tk = 3, figure
  for mj = [1 3] % 1:size(Ft,2)
    fvp = Ft.mjdg==mj&Ft.task==tk&~rt4s_out&Ft.qu1; 
    fvn = Ft.mjdg==mj&Ft.task==tk&~rt4s_out&~Ft.qu1; 
    lgc = {'Yes' 'No'};  % {'Correct' 'Incorrect'}; {'Yes' 'No'}
    if tk==2 
        if islr
            fvp = Ft.mjdg==mj&Ft.task==tk&~rt4s_out&Ft.tk2qu1lr;
            fvn = Ft.mjdg==mj&Ft.task==tk&~rt4s_out&~Ft.tk2qu1lr; 
            lgc = {'CCW' 'CW'};
        else
            fvp = Ft.mjdg==mj&Ft.task==tk&~rt4s_out&Ft.right1; 
            fvn = Ft.mjdg==mj&Ft.task==tk&~rt4s_out&Ft.wrong1; 
            lgc = {'Correct' 'Incorrect'}; 
        end
    elseif tk==3
        lgc = {'Correct' 'Incorrect'}; 
    end

    for sp = [1 3]
        subplot(2, 2, sp+(mj-1)/2)
        axis([xmin(tk) xmax(tk) ymin(tk) ymax(tk)]), hold on, box on, grid on
        if sp == 1
            [m, s, n, g] = grpstats(Ft{fvp,'rti1'}, Ft{fvp,'ftr1'}, {'mean', 'sem', 'numel', 'gname'});
            eb(1)=errorbar(1:length(m), m, s);
            [m, s, n, g] = grpstats(Ft{fvn,'rti1'}, Ft{fvn,'ftr1'}, {'mean', 'sem', 'numel', 'gname'});
            eb(2)=errorbar(1:length(m), m, s);
            if     tk==1 && mj == 1, title('Detection')
            elseif tk==1 && mj == 3, title('Control')     
            elseif tk==2 && mj == 1, title({'Categorization' '(CR session)'})
            elseif tk==2 && mj == 3, title({'Categorization' '(FVR session)'})  
            elseif tk==3 && mj == 1, title({'Identification' '(CR session)'})
            elseif tk==3 && mj == 3, title({'Identification' '(FVR session)'})              
            end
        else
            [m, s, n, g] = grpstats(Ft{fvp,'rti2'}, Ft{fvp,'ftr1'}, {'mean', 'sem', 'numel', 'gname'});
            eb(1)=errorbar(1:length(m), m, s);
            [m, s, n, g] = grpstats(Ft{fvn,'rti2'}, Ft{fvn,'ftr1'}, {'mean', 'sem', 'numel', 'gname'});
            eb(2)=errorbar(1:length(m), m, s);
            xlabel('Signal strength')  
            if     mj == 1, title('Confidence judgment')
            elseif mj == 3, title('Visibility judgment')     
            end
        end
        if mj==1
            ylabel('Decision time (s)')
            legend(eb, lgc, 'location', 'best')
        elseif tk~=1, legend(eb, lgc, 'location', 'best')        
        end  
        ca = gca;
        set(ca.Children(1), 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '--')
        set(ca.Children(2), 'LineWidth', 1.5, 'Color', 'k')
        %set(ca, 'XTickLabel', {})
    end    
  end
end

subplot(221), boxplot(Ft{Ft.task==1 & ~rt4s_out, 'rti1'}, Ft{Ft.task==1 & ~rt4s_out, 'ftr1'})
subplot(222), boxplot(Ft{Ft.task==1 & ~rt4s_out, 'rti2'}, Ft{Ft.task==1 & ~rt4s_out, 'ftr1'})
subplot(223), boxplot(Ft{Ft.task==2 & ~rt4s_out, 'rti1'}, Ft{Ft.task==2 & ~rt4s_out, 'ftr1'})
subplot(224), boxplot(Ft{Ft.task==2 & ~rt4s_out, 'rti2'}, Ft{Ft.task==2 & ~rt4s_out, 'ftr1'})

% boxplots
figure
boxplot(Ft.ps_offset_st, {Ft.qu2 Ft.mjdg}, 'jitter', 0.1, 'notch', 'on', ...
    'datalim', [-inf inf], 'extrememode', 'compress', 'whisker', 1.5)
figure
[m, s, n, g] = grpstats(Ft.qu2, Ft.mjdg, {'mean', 'sem', 'numel', 'gname'})
errorbar(1:length(m), m, s)

% mjdg-task-ftr1 histograms
for m = [1 3]
  for i = 1:3
    % waterfall, meshc, surfc, contourf, contour3, surface
    %figure(1)
    %subplot(2,3,3/2*(m-1)+i), histogram2(Ft{Ft.mjdg==m & Ft.task==i, 'qu1'}, ...
    %    Ft{Ft.mjdg==m & Ft.task==i, 'ftr1'}) 
    
    figure(2)
    subplot(2,3,3/2*(m-1)+i), histogram2(Ft{Ft.mjdg==m & Ft.task==i, 'qu12'}, ...
        Ft{Ft.mjdg==m & Ft.task==i, 'ftr1'}, 'FaceAlpha', 0.7) 
    
    %figure(3)
    %subplot(2,3,3/2*(m-1)+i), histogram2(Ft{Ft.mjdg==m & Ft.task==i, 'qu2'}, ...
    %    Ft{Ft.mjdg==m & Ft.task==i, 'ftr1'}) 
  end
end

% correlations and partial correlations
subplot(131), surfc(corr(Ft{:,[1 2 4:6 8 9 11:28]}))
subplot(132), contourf(abs(corr(Ft{:,[1 2 4:6 8 9 11:28]})))
subplot(133), surfc(partialcorr(Ft{:,[1 2 4:6 8 9 11:28]}))

% scatterplots
figure
gscatter(Ft.rti2, Ft.ps_onset_qu2, Ft.cm, [], 'oxxo', 6)
figure
scatterhist(Ft.rti1, Ft.rti2)
figure
hist3([Ft.rti1 Ft.rti2], [40 40])
figure
gplotmatrix(Ft{:, {'ftr1' 'ps_onset_qu2'}}, Ft{:, {'ps_rt1' 'ps_onset_qu2'}}, ...
    Ft.cm, [], 'oxxo', 4, 'on', 'hist', {'ftr1' 'ps onsetqu2'}, {'ps rt1' 'ps onsetqu2'})

gplotmatrix(Ft{:, {'ftr1' 'rti1'}}, Ft{:, {'rti2' 'qu2'}}, Ft.cm)


% fig05bottom, fig08e: Between subject scatterplots da vs meta_da
tk = 1; 
if     tk==1, Lim=[-.5 4.5];  Tick=0:4;
elseif tk==2, Lim=[-.5 2.75]; Tick=0:.5:3; end
for i = [1 3]
  if tk==1
      figure
      G = G4f13(G4f13.mjdg==i & G4f13.task==tk & G4f13.ftr1~=1, :);
      gplotmatrix(G.da, G.meta_da, G.ftr1, ...
          repmat(linspace(0.9,0.2,11-tk*4),3,1)','.',20,'on',[],'d_a','meta-d_a')
  elseif tk==2
      G = G4f13(G4f13.task==tk & G4f13.ftr1~=1, :);
      gplotmatrix(G.da, G.meta_da, {G.ftr1 G.mjdg}, ...
          kron(linspace(0.9,0.2,11-tk*4),ones(3,2))','.*',[20 10],'on',[],'d_a','meta-d_a')
  end
  cf = gcf;
  set(cf.Children(1), 'Location', 'southeast')
  lgs = cf.Children(1).String; 
  if tk==2, lgs = {'L2 CR' 'L2 FVR' 'L3 CR' 'L3 FVR' 'L4 CR' 'L4 FVR'}; end
  set(cf.Children(1), 'String', lgs)
  set(cf.Children(2),'XGrid', 'on', 'YGrid', 'on')
  set(cf.Children(2), 'XLim', Lim, 'YLim', Lim, 'YTick', Tick)  
  refline(cf.Children(2), 1, 0)
  set(cf.Children(2).Children(1), 'Color', 'k', 'LineStyle', ':', 'HandleVisibility', 'off')    
end


% between-subjects CR-VR scatterplot of qu2 by ftr1 and task
T1 = Ft4(Ft4.task==1 & Ft4.qu1==1 & Ft4.ftr2~=categorical(2),:);
ts1 = grpstats(T1, {'mjdg' 'ftr1' 'sid'}, {'mean'}, 'DataVars', 'qu2');
U1 = unstack(ts1, {'mean_qu2'}, 'mjdg', 'GroupingVariables', {'ftr1' 'sid'}, 'NewDataVariableNames', {'CR' 'VR'});
gscatter(U1.CR, U1.VR, U1.ftr1, [], 'o', 5)

T2 = Ft4(Ft4.task~=1 & Ft4.ftr2~=categorical(2),:);
ts2 = grpstats(T2, {'mjdg' 'task' 'ftr1' 'sid' 'qu1'}, {'mean'}, 'DataVars', 'qu2');
U = unstack(ts2, {'mean_qu2'}, 'mjdg', 'GroupingVariables', {'task' 'ftr1' 'sid' 'qu1'}, 'NewDataVariableNames', {'CR' 'VR'});
%T2 = Ft4(Ft4.mjdg==1 & Ft4.ftr2~=categorical(2),:);
%ts2 = grpstats(T2, {'task' 'ftr1' 'sid'}, {'mean'}, 'DataVars', {'qu1' 'qu2'});

subplot(1,2,1), axis([1 4 1 4]), hold on    
gscatter(U{U.task==2,'CR'}, U{U.task==2,'VR'}, U{U.task==2,'qu1'}, [], 'ox', 8)
subplot(1,2,2), axis([1 4 1 4]), hold on    
gscatter(U{U.task==2,'CR'}, U{U.task==2,'VR'}, U{U.task==2,'ftr1'}, [], 'o', 6:14)
%gscatter(U{U.task==i,'CR'}, U{U.task==i,'VR'}, U{U.task==i,'sid'}, [], 'ooooox', 10)


% between-subjects CR-VR scatterplot of meta_da,t1ca,da by discr task
%G = G(G.task~=1 & G.ftr2~=2,:); cstr={'da' 'meta_da' 'c_1'};
G = G4f13(G4f13.task==1,:); cstr={'da' 'meta_da' 'c_1'};
for i=1:3 
  subplot(1,3,i)
  U = unstack(G, cstr(i), 'mjdg', 'GroupingVariables', {'task' 'ftr1' 'sid'}, 'NewDataVariableNames', {'CR' 'VR'});
  gscatter(U.CR, U.VR, U.ftr1, [], 'o', 4:12)
end
gscatter(U.CR, U.VR, U.ftr1, [], 'o', 4:12)


%between subjects guessing criterion
T3 = Ft4(any(Ft4.task==[1 3],2) & Ft4.ftr2~=categorical(2),:); 
ts3 = grpstats(T3, {'task' 'mjdg' 'qu1' 'qu2' 'sid'}, {'mean'}, 'DataVars', 'right1');
unstack(ts3, {'mean_right1'}, 'qu2', 'GroupingVariables', {'task' 'mjdg' 'qu1'}, 'AggregationFunction', @mean, 'NewDataVariableNames', {'R1' 'R2' 'R3' 'R4'})
unstack(ts3, {'mean_right1'}, 'qu2', 'GroupingVariables', {'task' 'mjdg' 'qu1'}, 'AggregationFunction', @std, 'NewDataVariableNames', {'R1' 'R2' 'R3' 'R4'})

T3 = Ft4(Ft4.task==2 & Ft4.ftr2~=categorical(2),:);
T3.lr = double(T3.ftr1_raw > 4);
ts3 = grpstats(T3, {'mjdg' 'lr' 'qu2' 'sid'}, {'mean'}, 'DataVars', 'qu1');
unstack(ts3, {'mean_qu1'}, 'qu2', 'GroupingVariables', {'mjdg' 'lr'}, 'AggregationFunction', @mean, 'NewDataVariableNames', {'R1' 'R2' 'R3' 'R4'})
unstack(ts3, {'mean_qu1'}, 'qu2', 'GroupingVariables', {'mjdg' 'lr'}, 'AggregationFunction', @std, 'NewDataVariableNames', {'R1' 'R2' 'R3' 'R4'})

% pooled guessing criterion
T3 = Ft4(any(Ft4.task==[1 3], 2) & Ft4.ftr2~=categorical(2),:); 
grpstats(T3, {'task' 'mjdg' 'qu1' 'qu2'}, {'mean' 'sem'}, 'DataVars', 'right1')
T3 = Ft4(Ft4.task==2 & Ft4.ftr2~=categorical(2),:); T3.lr = double(T3.ftr1_raw > 4);
grpstats(T3, {'mjdg' 'lr' 'qu2'}, {'mean' 'sem'}, 'DataVars', 'qu1')

