function [logdethr,q] = E3calib_det(pThreshold,viswon,noam,nt,spo,lang)
% E3_4: calibration of stimulus absolute detection(Y/N) threshold 
% Block of trials serving as calibration task
% Calibration method: Quest adaptive staircase
% Monotropical: 90


KbName('UnifyKeyNames');
qu1keys{1}='2'; % Y
qu1keys{2}='3'; % N
RestrictKeysForKbCheck([KbName(qu1keys{1}) KbName(qu1keys{2}) KbName('space')]);
oldvdl=Screen('Preference','VisualDebugLevel',1);
show_results=0;
LoadAcerAF715GammaTable=1;
debugmode=0;

if viswon  %this switch does not work
    PsychVideoSwitcher('SwitchMode',0,viswon,0); %high precision luminance mode
end


%%%%%%%%%%%%%%%%%
% initialization
%%%%%%%%%%%%%%%%%

% Defining parameters: layout
if viswon, par.gsb=14;         % grayscale shades bits
else       par.gsb=8;
end
par.mgv=2^par.gsb-1;    % maximum gun value: number of shades of gray 
par.btrr=126.3;
par.fplu=0.5; par.fplucl=par.fplu*par.mgv;   % fixation point luminance
par.texlu=0.6;                     % text luminance
% Defining parameters: stimuli
par.lambda=120;            % gabor sine grating spatial frequency 
par.bglu=0.4; par.bglucl=par.bglu*par.mgv;   % stimuli background luminance
par.noam=noam;  par.noamcl=par.noam*par.mgv;  % noise amplitude 
%The mean and variance are always specified as if the image were of class double in the range [0 1]. 
%If the input image is of class uint8 or uint16, the imnoise function converts the image to double, 
%adds noise according to the specified type, and then converts the noisy image back to the same class as the input
% Defining paremeters: trial sequence
par.nt=nt;                    % total number of trials
par.nq=1;                         % number of questions per trial   
par.presti=1;           % prestimulus interval: baseline for PS
par.presi0=.3;par.presi1=.3;         % prestimulus intervals
par.sdi=.3;                   % stimulus display time in ms
par.posti=0;                     % poststimulus interval
par.respwi=inf;                   % response window
par.postresp2=1;                  % window for measuring PS after last question
% ITI=presi0+presi1*rand+stimdisp+posti+respwi=300+400+200:300+300+400+200+1200=900:2400


% Definition of PF intensity: lwc=log10(talu/bglu)
%logarithm of Weber contrast of Gabor peak amplitude 
%talu=[1:128], lwc[log10(0.78%:100%)=-2.1:0]
par.stluGuess=0.001; %log10(0.001)=-3

% Provide our prior knowledge to QuestCreate, and receive the data struct q
%global pThreshold logdethr logdetsd q
tGuess=log10(par.stluGuess/par.bglu);   % [-5:1] 
tGuessSd=3;
%pThreshold=0.75;                   % threshold criterion expressed as P(response)==1           
beta=3.5;delta=0.01;
gamma=0.5;                         % because it is a Y/N task, gamma=0.5
grain=0.01;range=6;                % intensity as logarithm of contrast
par.pvlurgmin=0;          % minimum gabor peak-valley amplitude allowed 
par.pvlurgmax=0.4;          % maximum gabor peak-valley amplitude allowed 
par.brlosd=0.08;
q=QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range);
q.normalizePdf=1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials

% Preparing FP frame
[~,par.l,c,par.csrad,par.sigma]=E3gabor_noise(par.bglu,par.noam,par.lambda,45,5);       % extracting image of size l
EM64=par.bglu*ones(par.l);         % empty matrix 
FM64=EM64;   %FM64(c,c)=par.fplu;% empty matrix with simple fixation point 
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


% Image sequence presentation and control flow statements

% Instructions text
if ~isnan(str2double(qu1keys{1}))
    padkeys{1}=['Pad' qu1keys{1}]; padkeys{2}=['Pad' qu1keys{2}];
end
if strcmp(lang,'EN')
    instr{1}='Calibration task';
    instr{2}=' '; %'This task will measure your detection threshold.';
    instr{3}='Vertical gratings will be presented sometimes at the center of the screen.';
    instr{4}='By selecting the corresponding keys, you will have to report:';
    instr{5}='Detection';
    instr{6}=['[' padkeys{1} ']: Seen(Y)   [' padkeys{2} ']: Not-seen(N)']; % instr{6}='[Y]: detection   [N]: non-detection';
    instr{7}='IMPORTANT: answer as fast as possible, but give priority to response accuracy';
    instr{8}='Press spacebar when ready';
    instr{9}='Press any key to continue';
    instr{10}='{Y N}';
    switch OSName
        case 'Linux', fns='-misc-fixed-medium-r-normal--18-120-100-100-c-90-iso8859-15';
        case 'Windows', fns='Helvetica';
        otherwise, fns='Helvetica';
    end
elseif strcmp(lang,'JP')
    instr{1}='校正課題';
    instr{2}=' '; %この課題では、あなたの検知能力を計ります
    instr{3}='画面の中央に縦の縞模様が表示されます。';
    instr{4}='対応するキーを押すことで、次の模様の特徴を報告してください。';
    instr{5}='検知';
    instr{6}=['[' padkeys{1} ']: 見えた(Y)   [' padkeys{2} ']: 見えなかった(N)']; 
    instr{7}='できるだけ早く答えて下さい。ただし、早さより正確さを優先してください。';
    instr{8}='準備ができたら、スペースキーを押して始めてください。';
    instr{9}='次に進んでもよければ任意のキーを押してください';
    instr{10}='{Y N}';
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
    par.filename=['E3calib_det_',dateISO8601];
    his.tint=[];
    his.kbi=[];
    his.det=[];
    %logfid=fopen(['datalog_',dateISO8601,'.txt'],'w');
    %fprintf(logfid,'Start time: %s\n',dateISO8601);
    %resfid=fopen(['datares_',dateISO8601,'.txt'],'w');
    his.timelog=[];

    % Stacking frames in offscreen buffers 
    wip(7)=Screen(0,'OpenOffscreenWindow',par.gcl);
    Screen(wip(7),'TextColor',par.texlu*par.mgv*[1 1 1]);
    Screen(wip(7),'TextFont',fns);
    Screen(wip(7),'TextSize',40); bourec=Screen('TextBounds', wip(7),instr{1});
    Screen(wip(7),'DrawText',instr{1},(par.scrsiz(3)-bourec(3))/2,0);
    Screen(wip(7),'TextSize',20); bourec=Screen('TextBounds', wip(7),instr{2});
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
    %tstring='Calibration task\n\n\nThis task will measure your detection threshold.\nNoisy vertical gratings will be presented at the center of the screen.\nBy selecting the corresponding keys, you will have to report:\n\nDetection\n[Y]: detection   [N]: non-detection\n\n\n\nPress spacebar when ready';
    %DrawFormattedText(wip(7), tstring,'center','center',par.texlu*par.mgv*[1 1 1]);

    wip(2)=Screen(0,'OpenOffscreenWindow',par.gcl);
    Screen(wip(2),'TextColor',par.texlu*par.mgv*[1 1 1]);
    Screen(wip(2),'TextFont',fns); 
    bourec=Screen('TextBounds', wip(2),instr{9});
    Screen(wip(2),'DrawText',instr{9},(par.scrsiz(3)-bourec(3))/2,.9*par.scrsiz(4));

    %wip(3)=Screen(0,'OpenOffscreenWindow',par.gcl);Screen(wip(3),'PutImage',FM);
    wip(3)=Screen(wip(1),'MakeTexture',FM);

    %wip(5)=Screen(0,'OpenOffscreenWindow',par.gcl);Screen(wip(5),'PutImage',EM);        % blank 
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

    Screen(wip(1),'DrawTexture',wip(7));
    [~, tt]=Screen('Flip',wip(1));

    %while GetChar ~=' ', FlushEvents; WaitSecs(0.1); end, t0=GetSecs; 
    [t0,KeyCode]=KbWait([],2);
    while ~KeyCode(KbName('space'))==1, [t0,KeyCode]=KbWait([],2); end
    if isobject(spo) && strcmp(spo.Status,'open')
        spstr='t0';
        fprintf(spo,'ET_REM %s\n',spstr); 
    end
    
    erco=0;
    vatr=0; % valid trials
    
    while vatr < par.nt || QuestSd(q) > par.brlosd
        % Interleave half of trials with no stimulus, lest the subject assigns
        % high probability to stimulus presence in the main task
        if rand > 0.7  % stiabs=round(rand);
            stiabs=1;
        else
            stiabs=0;
        end

        tsug=QuestQuantile(q);	%recommended by Pelli (1987), and still our favorite

        %we are free to test any intensity we like, not necessarily what Quest suggested
        tint=par.bglu*10^tsug;
        tint=min(par.pvlurgmax,max(par.pvlurgmin,tint)); %restrict range and grain of contrasts
        tint=round(tint*par.mgv)/par.mgv; 
        % issue error if the display grayscale resolution is surpassed 3 times
        if tint < 1/par.mgv
            erco=erco+1;
            tint=1/par.mgv; % avoid overflowing the pdf
            if erco > 2
                error('KuvikException:ExceededDeviceResolution',...
                    'Monitor grayscale resolution is too coarse. Calibration aborted');
            end
        end
        tsug_r=log10(tint/par.bglu);
        % stimulus absent as control
        if stiabs
            tint=0;
            tsug_r=-inf;
        end

        %draw gabor into buffer, all parameters are passed on the fly
        I64=E3gabor_noise(par.bglu,par.noam,par.lambda,90,tint); 
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

        %wip(4)=Screen(0,'OpenOffscreenWindow',par.gcl);Screen(wip(4),'PutImage',I);
        wip(4)=Screen(wip(1),'MakeTexture',I);

        %{
        Screen('CopyWindow',wip(2),wip(1)); % pausing between trials
        pause_t=Screen('Flip',wip(1));                  
        [s,p]=KbWait([],1,10+pause_t);                                  
        if p==KbName('ESCAPE'), sca, return, end
        clear s p
        %}
        
        % Baseline window
        WaitSecs(par.presti);

        % presenting stimulus (padded prelaterally with random intervals)
        Screen(wip(1),'DrawTexture',wip(3)); [~, onset_fp]=Screen('Flip',wip(1));
        %fprintf(logfid,'%d: stimulus %d on\n',onset_fp,i);
        WaitSecs(par.presi0+par.presi1*rand);
        if debugmode
            Screen('CopyWindow',wip(4),wip(1));[~, onset_st]=Screen('Flip',wip(1));
        else
            Screen(wip(1),'DrawTexture',wip(4));[~, onset_st]=Screen('Flip',wip(1));
        end
        %fprintf(logfid,'%d: stimulus %d off\n',offset_fp,i);
        if isobject(spo) && strcmp(spo.Status,'open')
            spstr=['onset_st'];
            fprintf(spo,'ET_REM %s\n',spstr); 
        end
        WaitSecs(par.sdi);
        Screen(wip(1),'DrawTexture',wip(5)); [~, offset_st]=Screen('Flip',wip(1));
        WaitSecs(par.posti);

        % recording subject feedback                
        Screen(wip(1),'DrawTexture',wip(6)); [~, resp_t]=Screen('Flip',wip(1));
        [rt,KeyCode]=KbWait([],2,par.respwi+resp_t);
        if ~isempty(find(KeyCode,1)) && length(find(KeyCode))==1
            kbi=find(KeyCode);
            if find(KeyCode)==KbName(qu1keys{1})
                whichdet=1;
            elseif find(KeyCode)==KbName(qu1keys{2})
                whichdet=0;
            else
                whichdet=nan;
            end
        else       
            kbi=nan;
            whichdet=nan;
        end
        KbReleaseWait; 
        
        % Wait after FOJ to allow PS to stabilize
        Screen(wip(1),'DrawTexture',wip(5)); Screen('Flip',wip(1));
        WaitSecs(par.postresp2);
        
        his.timelog=[his.timelog [tt-t0;t0-t0;onset_fp-t0;onset_st-t0;...
            offset_st-t0;rt-t0]];

        if isempty(whichdet) || isnan(whichdet) || stiabs
            continue
        else % update the pdf
            % Add the new datum (test intensity, observer response) to database
            q=QuestUpdate(q,tsug_r,whichdet); 
            % writing in results structure
            his.tint=[his.tint [tsug;tsug_r;tint]];
            his.kbi=[his.kbi kbi];
            his.det=[his.det whichdet];
            %fprintf(resfid,'%d\n',c_det) 
            whichdet=[];
            vatr=vatr+1;
        end    
    end

    Screen('CloseAll');
    Screen('Preference','VisualDebugLevel',oldvdl);
    Screen('LoadNormalizedGammaTable',0,oldgt); % default: oldgt=repmat((0:1/255:1)',1,3) 
    RestrictKeysForKbCheck([]);
    ShowCursor;
    Priority(0);

    % Ask Quest for the final estimate of threshold
    logdethr=QuestMean(q); % Recommended by Pelli(1989) and King-Smith etal.(1994). Still our favorite.
    logdetsd=QuestSd(q);

    % Save data
    par.timesegment(2,:)=datevec(now);
    if isobject(spo) && strcmp(spo.Status,'open')
        spstr='t_end';
        fprintf(spo,'ET_REM %s\n',spstr); 
    end
    par.logdethr=logdethr; par.logdetsd=logdetsd;
    save(par.filename,'par','his','q');
    %fclose(logfid,resfid)


    if show_results
        % Printing results
        fprintf('The threshold estimate in log10 of contrast (mean+-sd): %.2f ?????????ｽ?????ｽ???????ｽ?????ｽ????????ｽ?????ｽ???????ｽ?????ｽ} %.2f\n',logdethr,logdetsd);
        pcdethr=100*10^logdethr;
        pcsdleft=100*(10^(logdethr)-10^(logdethr-logdetsd)); 
        pcsdright=100*(10^(logdethr+logdetsd)-10^(logdethr));  
        fprintf(1,'in %%contrast (mean+[-leftsd,+rightsd): %.2f+[-%.2f,+%.2f]\n',pcdethr,pcsdleft,pcsdright);
        cldethr=(par.bglu*par.mgv/100)*pcdethr;
        clsdleft=(par.bglu*par.mgv/100)*pcsdleft; 
        clsdright=(par.bglu*par.mgv/100)*pcsdright;  
        fprintf('in target amplitude RGB %d-bit luminance (mean+[-leftsd,+rightsd): %.2f+[-%.2f\n,+%.2f]\n',par.gsb,cldethr,clsdleft,clsdright);

        % Plotting results
        %sorted list of intensities and response frequencies. t=QuestTrials(q,0.1)
        t=QuestTrials(q); 
        fprintf(' intensity     p fit         p    trials\n');
        disp([t.intensity; QuestP(q,t.intensity-logdethr);(t.responses(2,:)./sum(t.responses)); sum(t.responses)]');

        % (possibly unnormalized) probability density of candidate thresholds
        x=log10(2^(1-par.gsb)):grain:log10(1);
        figure, plot(x,QuestPdf(q,x)), title('(possibly unnormalized) pdf of candidate thresholds')  

        % psychometric function: PF=QuestP(q,x) 
        % probability of a correct (or yes) response at intensity x, assuming threshold is at x=0
        x=-1:0.01:1;
        figure, plot(x,QuestP(q,x)),title(sprintf('PF at intensity x, assuming x=0 is at threshold=%.2f',logdethr))
    end

catch ME
    
    Screen('CloseAll');
    Screen('Preference','VisualDebugLevel',oldvdl);
    Screen('LoadNormalizedGammaTable',0,oldgt); %oldgt is typically repmat((0:1/255:1)',1,3)
    RestrictKeysForKbCheck([]);
    ShowCursor;
    Priority(0);
    
    rethrow(ME);
end
