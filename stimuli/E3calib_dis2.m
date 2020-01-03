function [logdisth,logdissd,q]=E3calib_dis2(pThreshold,viswon,noam,nt,spo,lang)
% E3_4: calibration of stimulus orientation(4AFC) discrimination threshold 
% Block of trials serving as calibration task
% Calibration method: Quest adaptive staircase
% Quadratropical: 0,45,90,135


KbName('UnifyKeyNames');
RestrictKeysForKbCheck([KbName('o') KbName('9(') KbName('8*') KbName('7&')...
    KbName('6') KbName('9') KbName('8') KbName('7') KbName('space')]);
oldvdl=Screen('Preference','VisualDebugLevel',1);
show_results=0;
LoadAcerAF715GammaTable=1;
debugmode=0;

if viswon
    PsychVideoSwitcher('SwitchMode',0,viswon,0); %high precision luminance mode
end


%%%%%%%%%%%%%%%%%
% initialization
%%%%%%%%%%%%%%%%%

% Defining paremeters: layout
if viswon, par.gsb=14;         % grayscale shades bits
else       par.gsb=8;
end
par.mgv=2^par.gsb-1;    % maximum gun value: number of shades of gray 
par.btrr=126.3;
par.fplu=0.5; par.fplucl=par.fplu*par.mgv;   % fixation point luminance
par.texlu=0.6;                     % text luminance
% Defining paremeters: stimuli
par.lambda=120;            % gabor sine grating spatial frequency 
par.bglu=0.4; par.bglucl=par.bglu*par.mgv;    % stimuli background luminance
par.noam=noam; par.noamcl=par.noam*par.mgv;    % noise amplitude
% The mean and variance are always specified as if the image were of class double in the range [0 1]. 
% If the input image is of class uint8 or uint16, the imnoise function converts the image to double, 
% adds noise according to the specified type, and then converts the noisy image back to the same class as the input
% Defining paremeters: trial sequence
par.no=4;
par.nt=nt;                    % total number of trials
par.nq=1;                         % number of questions per trial
par.presti=1;           % prestimulus interval: baseline for PS
par.presi0=.3;par.presi1=.3;         % prestimulus intervals
par.sdi=.3;                   % stimulus display time in s
par.posti=0;                     % poststimulus interval
par.respwi=inf;                   % response window
par.postresp2=1;                  % window for measuring PS after last question
%ITI=presi0+presi1*rand+sdi+posti+respwi=300+400+200:300+300+400+200+1200=900:2400


% Definition of PF intensity: lwc=log10(talu/par.bglu)
% logarithm of Weber contrast of Gabor peak amplitude 
% talu=[1:128], lwc[log10(0.78%:100%)=-2.1:0]
par.stluGuess=0.004; %log10(.1)=-2

% Provide our prior knowledge to QuestCreate, and receive the data struct q
tGuess=log10(par.stluGuess/par.bglu);   % [-2.1:0] 
tGuessSd=3;
%pThreshold2=0.625;                   % threshold criterion expressed as P(response)==1           
beta=3.5;delta=0.01;
gamma=0.25;                         % because it is a 4AFC task, gamma=0.25 (assuming a perceptual(vs subjective awareness) task)
grain=0.01;range=6;                % intensity as logarithm of contrast
par.pvlurgmin=0;          % minimum gabor peak-valley amplitude allowed 
par.pvlurgmax=0.4;          % maximum gabor peak-valley amplitude allowed
par.brlosd=0.08;
q=QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range);
q.normalizePdf=1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials

% Preparing FP frame
[~,par.l,c,par.csrad,par.sigma]=E3gabor_noise(par.bglu,par.noam,par.lambda,45,5);       % extracting image of size l
EM64=par.bglu*ones(par.l);         % empty matrix 
FM64=EM64;
FM64(c-3:c+3,c)=par.fplu;FM64(c,c-3:c+3)=par.fplu;     % crosshair FP 
if viswon
    EM14=PsychVideoSwitcher('MapLuminanceToRGB',EM64,par.btrr);
    FM14=PsychVideoSwitcher('MapLuminanceToRGB',FM64,par.btrr);
    EM=EM14; 
    FM=FM14;
else
    EM8=uint8(EM64*par.mgv);
    FM8=uint8(FM64*par.mgv);
    EM=EM8;
    FM=FM8;
end


% Shuffling algorithm: creating a random sequence of par.no categories
%o_rpseq=randperm(par.nt);             % orientation values sequence random permutation              
%o_rpseq=mod(o_rpseq,par.no)+1;      % shuffle instead of randomize (j=ceil(no*rand)) 
%for i=1:par.no
%    o_rpseq(o_rpseq==i)=45*(i-1); % 000,045,090,135
%end
%par.o_rpseq=o_rpseq;


% Image sequence presentation and control flow statements

% Instructions text
if strcmp(lang,'EN')
    instr{1}='Calibration task';
    instr{2}=' '; %'This task will measure your orientation discrimination threshold.';
    instr{3}='Gratings with different orientations will be presented at the center of the screen.';
    instr{4}='By selecting the corresponding keys, you will have to report:';
    instr{5}='Orientation discrimination (000,045,090,135)';
    instr{6}='[Pad6]: 000 degrees   [Pad9]: 045 degrees   [Pad8]: 090 degrees   [Pad7]: 135 degrees';
    instr{7}='IMPORTANT: answer as fast as possible, but give priority to response accuracy';
    instr{8}='Press spacebar when ready';
    instr{9}='Press any key to continue';
    instr{10}='{000 045 090 135}';
    switch OSName
        case 'Linux', fns='-misc-fixed-medium-r-normal--18-120-100-100-c-90-iso8859-15';
        case 'Windows', fns='Helvetica';
        otherwise, fns='Helvetica';
    end
elseif strcmp(lang,'JP')
    instr{1}='校正課題';
    instr{2}=' ';
    instr{3}='異なる傾きの縞模様が画面の中央に表示されます。';
    instr{4}='対応するキーを押すことで、次の模様の特徴を報告してください。';
    instr{5}='傾き識別 (000,045,090,135)';
    instr{6}='[Pad6]: 000 度   [Pad9]: 045 度   [Pad8]: 090 度   [Pad7]: 135 度';
    instr{7}='できるだけ早く答えて下さい。ただし、早さより正確さを優先してください。';
    instr{8}='準備ができたら、スペースキーを押して始めてください。';
    instr{9}='次に進んでもよければ任意のキーを押してください';
    instr{10}='{000 045 090 135}';
    instr=cellfun(@double,instr,'uniformOutput',false);
    fns='TakaoExGothic';
    if IsLinux, fns='-:lang=ja'; end
end
par.instr=instr;

% Getting parameters
[par.dispsiz(1) par.dispsiz(2)]=Screen('DisplaySize',0);
par.scrnum=Screen('Screens'); % 0 is the main screen
par.scrfmr=Screen('FrameRate',0);
par.scrdep=get(0,'ScreenDepth'); % up to 32; %pixdeps=Screen('PixelSizes',0);
oldgt=Screen('ReadNormalizedGammaTable',0);
if LoadAcerAF715GammaTable
    load AcerAF715; graygammatable=gammaTable1;
    Screen('LoadNormalizedGammaTable',0, graygammatable*[1 1 1]);
end
[par.gammatable,par.dacbits,par.reallutsize]=Screen('ReadNormalizedGammaTable',0);
par.wcl=WhiteIndex(0); par.bcl=BlackIndex(0); %CLUT index at the current screen depth
par.gcl=par.bglu*(par.wcl+par.bcl);
par.scrppi=get(0,'ScreenPixelsPerInch'); % 96
par.scrsiz=get(0,'ScreenSize'); %[par.scrsiz(3) par.scrsiz(4)]=Screen('WindowSize', 0);
switch debugmode
    case 0, winsiz=[];
    case 1, winsiz=[+10 +10 par.scrsiz(3)/2-10 par.scrsiz(4)/2-10];
    case 2, PsychDebugWindowConfiguration([],.9), winsiz=[];
end
par.timesegment(1,:)=datevec(now);


try 
    
    % Initializing data structures
    wip=zeros(1,par.nq+7);
    wip(1)=Screen(0,'OpenWindow',[],winsiz);%w(1)=0;
    [par.monitorFlipInterval par.nrValidSamples par.stddev]=Screen('GetFlipInterval',wip(1));
    Screen(wip(1),'FillRect',par.gcl);
    KbName('UnifyKeyNames'); 
    dateISO8601=datestr(now,30);
    par.filename=['E3calib_dis2_',dateISO8601];
    his.tint=[];
    his.tazim=[];
    his.kbi=[];
    his.dis=[];
    his.timelog=[];

    % Stacking frames in buffers
    wip(7)=Screen(0,'OpenOffscreenWindow',par.gcl);
    Screen(wip(7),'TextColor',par.texlu*par.mgv*[1 1 1]);
    Screen(wip(7),'TextFont',fns);
    Screen(wip(7),'TextSize',40);bourec= Screen('TextBounds', wip(7),instr{1});
    Screen(wip(7),'DrawText',instr{1},(par.scrsiz(3)-bourec(3))/2,0);
    Screen(wip(7),'TextSize',20);bourec= Screen('TextBounds', wip(7),instr{2});
    Screen(wip(7),'DrawText',instr{2},(par.scrsiz(3)-bourec(3))/2,.2*par.scrsiz(4));
    bourec= Screen('TextBounds', wip(7),instr{3});
    Screen(wip(7),'DrawText',instr{3},(par.scrsiz(3)-bourec(3))/2,.25*par.scrsiz(4));
    bourec= Screen('TextBounds', wip(7),instr{4});
    Screen(wip(7),'DrawText',instr{4},(par.scrsiz(3)-bourec(3))/2,.35*par.scrsiz(4));
    Screen(wip(7),'TextSize',22);bourec= Screen('TextBounds', wip(7),instr{5});
    Screen(wip(7),'DrawText',instr{5},(par.scrsiz(3)-bourec(3))/2,.45*par.scrsiz(4));
    Screen(wip(7),'TextSize',20);bourec= Screen('TextBounds', wip(7),instr{6});
    Screen(wip(7),'DrawText',instr{6},(par.scrsiz(3)-bourec(3))/2,.5*par.scrsiz(4));
    bourec= Screen('TextBounds', wip(7),instr{7});
    Screen(wip(7),'DrawText',instr{7},(par.scrsiz(3)-bourec(3))/2,.7*par.scrsiz(4));
    bourec= Screen('TextBounds', wip(7),instr{8});
    Screen(wip(7),'DrawText',instr{8},(par.scrsiz(3)-bourec(3))/2,.9*par.scrsiz(4));

    %wip(2)=Screen(0,'OpenOffscreenWindow',par.gcl);
    %Screen(wip(2),'TextColor',par.texlu*par.mgv*[1 1 1]);
    %Screen(wip(2),'TextFont',fns); bourec= Screen('TextBounds', wip(2),'Press any key to continue');
    %Screen(wip(2),'DrawText','Press any key to continue',(par.scrsiz(3)-bourec(3))/2,.9*par.scrsiz(4));

    wip(3)=Screen(wip(1),'MakeTexture',FM);

    wip(5)=Screen(wip(1),'MakeTexture',EM);

    wip(6)=Screen(0,'OpenOffscreenWindow',par.gcl);
    Screen(wip(6),'TextColor',par.texlu*par.mgv*[1 1 1]);
    Screen(wip(6),'TextFont',fns);
    bourec= Screen('TextBounds', wip(6),instr{10});
    Screen(wip(6),'DrawText',instr{10},(par.scrsiz(3)-bourec(3))/2,.9*par.scrsiz(4));


%%%%%%%%%%%%%%%%%
% task
%%%%%%%%%%%%%%%%%
    
    HideCursor;
    par.prilev=MaxPriority(wip(1));
    Priority(par.prilev);

    Screen(wip(1),'DrawTexture',wip(7)); %Screen(wip(1),'DrawTexture',wip(7));
    [~, tt]=Screen('Flip',wip(1));

    [t0,KeyCode]=KbWait([],2);
    while ~KeyCode(KbName('space'))==1, [t0,KeyCode]=KbWait([],2); end
    if isobject(spo) && strcmp(spo.Status,'open')
        spstr='t0';
        fprintf(spo,'ET_REM %s\n',spstr); 
    end
    
    erco=0;
    vatr=0;
    i=1;

    while vatr < par.nt || QuestSd(q) > par.brlosd
        tsug=QuestQuantile(q);	% Recommended by Pelli (1987), and still our favorite

        % We are free to test any intensity we like, not necessarily what Quest suggested
        tint=par.bglu*10^tsug;
        tint=min(par.pvlurgmax,max(par.pvlurgmin,tint)); % Restrict to range and grain of log contrasts that our equipment can produce
        tint=round(tint*par.mgv)/par.mgv; 
        % issue warning if the display grayscale resolution is surpassed
        if tint < 1/par.mgv
           erco=erco+1;
           tint=1/par.mgv; % avoid overflowing the pdf
           if erco > 2
               error('KuvikException:ExceededDeviceResolution',...
                   'Monitor grayscale resolution is too coarse. Calibration aborted');
           end
        end
        tsug_r=log10(tint/par.bglu);
        
        % randomizing orientation        
        o_rpseq(i)=45*(randi(4)-1);  % 000,045,090,135
        % computing stimulus
        I64=E3gabor_noise(par.bglu,par.noam,par.lambda,o_rpseq(i),tint); 
        if ~isempty(find(I64>1,1)) || ~isempty(find(I64<0,1))
            warning('KuvikException:ExceededRange',...
                'Out of luminance range. Rectifying noise distribution.');
            I64=min(1,max(0,I64));    % rectify the imagtrix: a ptbug creates artifacts otherwise
        end
        if viswon
            I14=PsychVideoSwitcher('MapLuminanceToRGB',I64,par.btrr);        
            I=I14;
        else
            I8=uint8(I64*par.mgv);    
            I=I8;
        end

        wip(4)=Screen(wip(1),'MakeTexture',I);   
        
        % Baseline window
        WaitSecs(par.presti);

        % presenting stimulus (padded bilaterally with random intervals)
        Screen(wip(1),'DrawTexture',wip(3)); [~,onset_fp]=Screen('Flip',wip(1));
        WaitSecs(par.presi0+par.presi1*rand);
        Screen(wip(1),'DrawTexture',wip(4)); [~,onset_st]=Screen('Flip',wip(1));
        if isobject(spo) && strcmp(spo.Status,'open')
            spstr=['onset_st'];
            fprintf(spo,'ET_REM %s\n',spstr); 
        end
        WaitSecs(par.sdi);
        Screen(wip(1),'DrawTexture',wip(5)); [~,offset_st]=Screen('Flip',wip(1));
        WaitSecs(par.posti);

        % recording subject feedback                
        Screen(wip(1),'DrawTexture',wip(6)); [~, resp_t]=Screen('Flip',wip(1));
        [rt,KeyCode]=KbWait([],2,par.respwi+resp_t);
        if ~isempty(find(KeyCode,1)) && length(find(KeyCode))==1
            kbi=find(KeyCode);
            switch kbi
                case {KbName('6'),KbName('9'),KbName('8'),KbName('7'),...
                        KbName('o'),KbName('9('),KbName('8*'),KbName('7&')}
                    switch kbi 
                        case {KbName('6'),KbName('o')}, whichdis=0;  % 000
                        case {KbName('9'),KbName('9(')}, whichdis=45;  % 045
                        case {KbName('8'),KbName('8*')}, whichdis=90;  % 090
                        case {KbName('7'),KbName('7&')}, whichdis=135;  % 135
                    end
                    whichdis=double(whichdis==o_rpseq(i));
                otherwise, whichdis=nan;
            end
        else
            kbi=nan;
            whichdis=nan;
        end
        KbReleaseWait; 

        % Wait between questions to allow PS to stabilize
        Screen(wip(1),'DrawTexture',wip(5)); Screen('Flip',wip(1));
        WaitSecs(par.postresp2);
        
        his.timelog=[his.timelog [tt-t0;t0-t0;onset_fp-t0;onset_st-t0;...
            offset_st-t0;rt-t0]];

        if isempty(whichdis) || isnan(whichdis)
            continue       
        else % update the pdf
            q=QuestUpdate(q,tsug_r,whichdis); % Add the new datum (actual test intensity and observer response) to the database
            % writing in results structure
            his.tazim=[his.tazim o_rpseq(i)];
            his.tint=[his.tint [tsug;tsug_r;tint]];       
            his.kbi=[his.kbi kbi];
            his.dis=[his.dis whichdis];
            whichdis=[];
            vatr=vatr+1;
        end    
        i=i+1;
    end

    Screen('CloseAll');
    Screen('Preference','VisualDebugLevel',oldvdl);
    Screen('LoadNormalizedGammaTable',0,oldgt);
    RestrictKeysForKbCheck([]);
    ShowCursor;
    Priority(0);

    % Ask Quest for the final estimate of threshold
    logdisth=QuestMean(q);		% Recommended by Pelli (1989) and King-Smith et al. (1994). Still our favorite.
    logdissd=QuestSd(q);

    % Save data
    par.timesegment(2,:)=datevec(now);
    if isobject(spo) && strcmp(spo.Status,'open')
        spstr='t_end';
        fprintf(spo,'ET_REM %s\n',spstr); 
    end
    par.logdisth=logdisth;par.logdissd=logdissd;
    save(par.filename,'par','his','q')
    
    
    if show_results
    end

catch ME
    
    Screen('CloseAll');
    Screen('Preference','VisualDebugLevel',oldvdl);
    Screen('LoadNormalizedGammaTable',0,oldgt);
    RestrictKeysForKbCheck([]);
    ShowCursor;
    Priority(0);
    
    rethrow(ME);
end

