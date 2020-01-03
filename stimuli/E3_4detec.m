 function E3_4detec(qu2mod,viswon,noam,logdethr,q,adbglu,adsfr,spo,lang)
% E3_4: Y/N detection task with CR or PAS
% Constant stimuli paradigm
% Monotropical: 90
%
% nq(Y/N,mcq)*ns     *nr= nq*nt
% 02         *(1*3*8)*30= 02*720  
% Randomize for nt=3*8*30=720=16*9*5=120*6 
% 6 blocks of ns*nr=24*5=120

 
KbName('UnifyKeyNames');
qu1keys{1}='2'; % Y
qu1keys{2}='3'; % N
RestrictKeysForKbCheck([KbName(qu1keys{1}) KbName(qu1keys{2}) ...
    KbName('1!') KbName('2@') KbName('3#') KbName('4$') ...
    KbName('space') KbName('c')]);
oldvdl=Screen('Preference','VisualDebugLevel',1);
show_results=0;
LoadAcerAF715GammaTable=1;
debugmode=0; 
pilotmode=0;

if viswon
    PsychVideoSwitcher('SwitchMode',0,viswon,0); %high precision luminance mode
end
if pilotmode % Pilot:4 shuffled sequences of 66 trials each: invoke uintrandseq 
    par.nr=12;
    par.ndetb=4;  
    par.ndett=72; 
    load E3_4detec_rps_pilot4b
else % 6 shuffled sequences of 120 trials each: invoke uintrandseq 
    par.nr=30; 
    par.ndetb=6; 
    par.ndett=120;
    load E3_4detec_rps
end
par.nrb=par.nr/par.ndetb;  %  nr per block
par.l_rps=l_rps;par.a_rps=a_rps;


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
par.texlu=0.6;                     %text luminance
% Defining paremeters: stimuli
par.logdethr=logdethr;
par.adsfr=adsfr; %64 adding extra target spatial freq factor
par.lambda=120;              %gabor sine grating spatial frequency 
par.adbglu=adbglu;%159 adding extra target average lum factor
par.adbglucl=par.adbglu*par.mgv;
par.bglu=0.4; par.bglucl=par.bglu*par.mgv; %stimulus control background luminance
par.noam=noam; par.noamcl=par.noam*par.mgv;    % noise amplitude
gamma=0.5; grain=0.01;
par.dethr=par.bglucl*10^logdethr;
ivq=find(QuestP(q,-1:grain:1)>=gamma+0.05,1,'first');
vq=QuestP(q,ivq*grain-1);  % virgintile=20-quantile, QuestP(q,[-1 0 1])~[.5001 .75 .995]; 
loglumvec=logdethr*ones(1,8)+vq*[-inf,-14/20,-9/20,-3/20,3/20,9/20,15/20,21/20];
par.lumvec=par.bglu*10.^loglumvec;  % target luminance values
par.lumveccl=par.lumvec*par.mgv;
% Defining paremeters: trial sequence
par.nl=length(par.lumvec); %#points in the pmf: target peak luminance values (contrast)  
par.na=3;               %control, additional average luminance series, and additional spatial frequency series  
par.ns=1*par.na*par.nl;                    % different stimuli types
par.nt=par.ndett*par.ndetb;               % total number of trials nt=nr*ns=par.ndett*par.ndetb
par.nq=2;                         % number of questions per trial   
par.presti=1;           % prestimulus interval: baseline for PS
par.presi0=.3; par.presi1=.3;         % prestimulus interval: randomizing onset time
par.sdi=.3;                   % stimulus display time in ms
par.posti=0;                     % poststimulus interval: allow immediate foj
par.intrespwi=1;             % gap to measure PS before soj
par.respwi1=inf;                   % response window for detection
par.respwi2=inf;                   % response window for CR
par.postresp2=1;                  % window for measuring PS after soj
%ITI=presi0+presi1*rand+stimdisp+posti+respwidet+respwicr=300+400+200:300+300+400+200+1200+3000=900:5400


% Stacking stimuli frames in buffer cells
[~,par.l,c,par.csrad,par.sigma]=E3gabor_noise(par.bglu,par.noam,par.lambda,45,5);% extracting image of size l
EM64=par.bglu*ones(par.l);         % empty matrix 
FM64=EM64;
FM64(c-3:c+3,c)=par.fplu;FM64(c,c-3:c+3)=par.fplu;     % crosshead FP 
BC64=cell(par.nl,par.na);             % cell holding nl*na images in buffer  
vadbglu=[1 adbglu/par.bglu 1];        % adding extra target average luminance set of trials
vadsfr=[1 1 adsfr/par.lambda];      % adding extra target spatial frequency set of trials
for i=1:par.nl 
    for j=1:par.na
        BC64{i,j}=E3gabor_noise(par.bglu*vadbglu(j),par.bglu,par.noam,par.lambda*vadsfr(j),90,par.lumvec(i)); 
    end
end

if viswon
    EM14=PsychVideoSwitcher('MapLuminanceToRGB',EM64,par.btrr);
    FM14=PsychVideoSwitcher('MapLuminanceToRGB',FM64,par.btrr);
    BC14=BC64;
    for i=1:par.nl 
        for j=1:par.na
            BC14{i,j}=PsychVideoSwitcher('MapLuminanceToRGB',BC64{i,j},par.btrr);
        end
    end
    EM=EM14; 
    FM=FM14; 
    BC=BC14;
else
    EM8=uint8(EM64*par.mgv);
    FM8=uint8(FM64*par.mgv);
    BC8=BC64;
    for i=1:par.nl 
        for j=1:par.na
            BC8{i,j}=uint8(BC64{i,j}*par.mgv);
        end
    end
    EM=EM8;
    FM=FM8; 
    BC=BC8;
end


% Image sequence presentation and control flow statements

% Instructions text
if ~isnan(str2double(qu1keys{1}))
    padkeys{1}=['Pad' qu1keys{1}]; padkeys{2}=['Pad' qu1keys{2}];
end
if strcmp(lang,'EN')
    instr{1}='Detection task';
    instr{2}=' '; %'This task will measure your perceptual sensitivity, and your subjective experience.';
    instr{3}='Vertical gratings will be presented at the center of the screen.';
    instr{4}='By selecting the corresponding keys, you will have to report:';
    instr{5}='Detection';
    if strcmp(qu2mod,'CR')
        instr{6}=['[' padkeys{1} ']: Seen(Y)   [' padkeys{2} ']: Not-seen(N)'];
        instr{7}='Confidence in your previous decision';
        instr{8}='[1]: not confident at all  [2]: a little confident  [3]: quite confident  [4]: absolutely confident';
    elseif strcmp(qu2mod,'PAS')
        instr{6}=['[' padkeys{1} ']: Seen(Y) (IMPORTANT: choose always detection as first answer here)'];
        instr{7}='Visibility of the stimulus';
        instr{8}='[1]: not seen  [2]: vaguely seen  [3]: almost clearly seen  [4]: absolutely clearly seen';
    else 
        error('KuvikException:InvalidArgument','Unrecognized function argument.');
    end
    instr{9}='IMPORTANT: try to use the {1 2 3 4} scale as exhaustively as possible';
    instr{10}='IMPORTANT: answer as fast as possible the first question after the stimulus, but give priority to response accuracy';
    instr{11}='You will receive your performance score at the end of the experiment';
    instr{12}='Press spacebar when ready';
    instr{13}='{Y N}';
    instr{14}='{1 2 3 4}';
    switch OSName
        case 'Linux', fns='-misc-fixed-medium-r-normal--18-120-100-100-c-90-iso8859-15';
        case 'Windows', fns='Helvetica';
        otherwise, fns='Helvetica';
    end
elseif strcmp(lang,'JP')
    instr{1}='検知課題';
    instr{2}=' '; %この課題では、あなたの知覚感度と主観的体験を計ります。
    instr{3}='画面の中央に縦の縞模様が表示されます。';
    instr{4}='対応するキーを押すことで、次の模様の特徴を報告してください。';
    instr{5}='検知課題';
    if strcmp(qu2mod,'CR')
        instr{6}=['[' padkeys{1} ']: 見えた(Y)   [' padkeys{2} ']: 見えなかった(N)'];
        instr{7}='直前の回答に対するあなたの確信度';
        instr{8}='[1]:全く自信がない  [2]:少し自信がある  [3]:かなり自信がある  [4]:完全に自信がある';
    elseif strcmp(qu2mod,'PAS')
        instr{6}=['[' padkeys{1} ']: 見えた(Y) 　重要：この課題では、常に「見えた」と報告してください'];
        instr{7}='縞模様の見え方';
        instr{8}='[1]:何も見えなかった  [2]:わずかに見えた  [3]:少し見えた  [4]:完全に見えた';
    else 
        error('KuvikException:InvalidArgument','Unrecognized function argument.');
    end
    instr{9}='重要：確信度を回答する際、課題に合わせてなるべく{1 2 3 4}スケールのすべての要素を使用して下さい。';
    instr{10}='重要：一番目の質問にできるだけ早く答えて下さい。ただし、早さより正確さを優先してください。';
    instr{11}='実験後、課題のスコアをお知らせします。';
    instr{12}='準備ができたら、スペースキーを押して始めてください。';
    instr{13}='{Y N}';
    instr{14}='{1 2 3 4}';
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
    % initialize KbCheck and variables to make sure they're properly allocated 
    % by Matlab - this to avoid time delays in the critical reaction time 
    % measurement part of the script
    [~,rt1,KeyCode1]=KbCheck;
    [~,rt2,KeyCode2]=KbCheck;

    % Initializing data structures
    wip=zeros(1,par.nq+7);
    wip(1)=Screen(0,'OpenWindow',[],winsiz);%w(1)=0;
    [par.monitorFlipInterval par.nrValidSamples par.stddev]=Screen('GetFlipInterval',wip(1));
    Screen(wip(1),'FillRect',par.gcl);
    KbName('UnifyKeyNames'); 
    dateISO8601=datestr(now,30);
    par.filename=['E3_4detec' qu2mod '_' dateISO8601];
    his.kbi=[];
    his.ftr=[];
    his.qu1=[];
    his.qu2=[];
    his.timelog=[];

    % Stacking frames in buffers
    wip(8)=Screen(0,'OpenOffscreenWindow',par.gcl);
    Screen(wip(8),'TextColor',par.texlu*par.mgv*[1 1 1]);
    Screen(wip(8),'TextFont',fns);
    Screen(wip(8),'TextSize',40); bourec=Screen('TextBounds',wip(8),instr{1});
    Screen(wip(8),'DrawText',instr{1},(par.scrsiz(3)-bourec(3))/2,0);
    Screen(wip(8),'TextSize',20); bourec= Screen('TextBounds', wip(8),instr{2});
    Screen(wip(8),'DrawText',instr{2},(par.scrsiz(3)-bourec(3))/2,.1*par.scrsiz(4));
    bourec=Screen('TextBounds',wip(8),instr{3});
    Screen(wip(8),'DrawText',instr{3},(par.scrsiz(3)-bourec(3))/2,.15*par.scrsiz(4));
    bourec=Screen('TextBounds',wip(8),instr{4});
    Screen(wip(8),'DrawText',instr{4},(par.scrsiz(3)-bourec(3))/2,.25*par.scrsiz(4));
    Screen(wip(8),'TextSize',22);bourec=Screen('TextBounds',wip(8),instr{5});
    Screen(wip(8),'DrawText',instr{5},(par.scrsiz(3)-bourec(3))/2,.35*par.scrsiz(4));
    Screen(wip(8),'TextSize',20);bourec=Screen('TextBounds',wip(8),instr{6});
    Screen(wip(8),'DrawText',instr{6},(par.scrsiz(3)-bourec(3))/2,.4*par.scrsiz(4));
    Screen(wip(8),'TextSize',22);bourec=Screen('TextBounds',wip(8),instr{7});
    Screen(wip(8),'DrawText',instr{7},(par.scrsiz(3)-bourec(3))/2,.5*par.scrsiz(4));
    Screen(wip(8),'TextSize',20);bourec=Screen('TextBounds',wip(8),instr{8});
    Screen(wip(8),'DrawText',instr{8},(par.scrsiz(3)-bourec(3))/2,.55*par.scrsiz(4));
    bourec= Screen('TextBounds',wip(8),instr{9});
    Screen(wip(8),'DrawText',instr{9},(par.scrsiz(3)-bourec(3))/2,.65*par.scrsiz(4));
    bourec= Screen('TextBounds',wip(8),instr{10});
    Screen(wip(8),'DrawText',instr{10},(par.scrsiz(3)-bourec(3))/2,.7*par.scrsiz(4));
    bourec= Screen('TextBounds',wip(8),instr{11});
    Screen(wip(8),'DrawText',instr{11},(par.scrsiz(3)-bourec(3))/2,.85*par.scrsiz(4));
    bourec=Screen('TextBounds',wip(8),instr{12});
    Screen(wip(8),'DrawText',instr{12},(par.scrsiz(3)-bourec(3))/2,.9*par.scrsiz(4));

    wip(3)=Screen(wip(1),'MakeTexture',FM);

    wip(5)=Screen(wip(1),'MakeTexture',EM);

    wip(6)=Screen(0,'OpenOffscreenWindow',par.gcl);
    Screen(wip(6),'TextColor',par.texlu*par.mgv*[1 1 1]);
    Screen(wip(6),'TextFont',fns);bourec=Screen('TextBounds',wip(6),instr{13});
    Screen(wip(6),'DrawText',instr{13},(par.scrsiz(3)-bourec(3))/2,.9*par.scrsiz(4));

    wip(7)=Screen(0,'OpenOffscreenWindow',par.gcl);
    Screen(wip(7),'TextColor',par.texlu*par.mgv*[1 1 1]);
    Screen(wip(7),'TextFont',fns);bourec=Screen('TextBounds',wip(7),instr{14});
    Screen(wip(7),'DrawText',instr{14},(par.scrsiz(3)-bourec(3))/2,.9*par.scrsiz(4));
    %Screen(wip(7),'DrawText','How confident are you about your decision? Select:{1 2 3 4}',(par.scrsiz(3)-bourec(3))/2,25);
    %Screen(wip(7),'DrawText','[K1]: not confident at all  [K2]: a little confident  [K3]: quite confident  [K4]: absolutely confident',(par.scrsiz(3)-bourec(3))/2,0);


%%%%%%%%%%%%%%%%%
% task
%%%%%%%%%%%%%%%%%
    
    HideCursor;
    par.prilev=MaxPriority(wip(1));
    Priority(par.prilev);

    Screen(wip(1),'DrawTexture',wip(8)); [~, tt]=Screen('Flip',wip(1));
    %Screen('CopyWindow',wip(8),wip(1)); [~, tt]=Screen('Flip',wip(1));
    
    [t0,KeyCode]=KbWait([],2);
    while ~KeyCode(KbName('space'))==1, [t0,KeyCode]=KbWait([],2); end
    if isobject(spo) && strcmp(spo.Status,'open')
        spstr='t0';
        fprintf(spo,'ET_REM %s\n',spstr); 
    end
    
    k=0;

    for j=1:par.ndetb
        WaitSecs(1);
        wip(9)=Screen(0,'OpenOffscreenWindow',par.gcl);
        Screen(wip(9),'TextColor',par.texlu*par.mgv*[1 1 1]);
        Screen(wip(9),'TextFont',fns);
        Screen(wip(9),'TextSize',30); bourec=Screen('TextBounds', wip(9),sprintf('Block %d',j));
        Screen(wip(9),'DrawText',sprintf('Block %d',j),(par.scrsiz(3)-bourec(3))/2,.3*par.scrsiz(4));
        Screen(wip(9),'TextSize',20);bourec=Screen('TextBounds', wip(9),instr{12});
        Screen(wip(9),'DrawText',instr{12},(par.scrsiz(3)-bourec(3))/2,.9*par.scrsiz(4));
        Screen(wip(1),'DrawTexture',wip(9));Screen('Flip',wip(1));
        [~,KeyCode]=KbWait([],2);
        while ~KeyCode(KbName('space'))==1 && ~KeyCode(KbName('c'))==1 
            [~,KeyCode]=KbWait([],2); 
        end
        Screen('Close',wip(9));
        if KeyCode(KbName('c'))==1
            iViewX2_SerialPort_Calib(wip(1),viswon);
        end
        
        for i=1:par.ndett
            % Putting stimulus in buffer 4
            wip(4)=Screen(wip(1),'MakeTexture',BC{l_rps(i,j),a_rps(i,j)});
            
            % Baseline window
            WaitSecs(par.presti);
            
            % presenting stimulus (padded prelaterally with random intervals)
            Screen(wip(1),'DrawTexture',wip(3)); [~, onset_fp]=Screen('Flip',wip(1));
            WaitSecs(par.presi0+par.presi1*rand);
            Screen(wip(1),'DrawTexture',wip(4)); [~, onset_st]=Screen('Flip',wip(1));   
            if isobject(spo) && strcmp(spo.Status,'open') 
                spstr='onset_st';%['E3_4detec',qu2mod,' onset_st-', num2str(l_rps(i,j)),'-',num2str(a_rps(i,j))];
                fprintf(spo,'ET_REM %s\n',spstr);  
            end
            % show stimulus until sdi elapses
            while (GetSecs-onset_st)<=par.sdi
                % poll for a response: subjects can respond before stimulus terminates
                if (KeyCode1(KbName(qu1keys{1}))==1 || KeyCode1(KbName(qu1keys{2}))==1)
                    break;
                end
                [~,rt1,KeyCode1]=KbCheck;
                % Wait 1 ms before checking the keyboard again to prevent overload 
                % of the machine at elevated Priority():
                WaitSecs(0.001);
            end
            WaitSecs('UntilTime',onset_st+par.sdi); % not necessary
            Screen(wip(1),'DrawTexture',wip(5)); [~, offset_st]=Screen('Flip',wip(1));
            WaitSecs(par.posti);

            % detection question: loop until valid key is pressed or reswi terminates
            Screen(wip(1),'DrawTexture',wip(6)); [~, onset_qu1]=Screen('Flip',wip(1));
            while (KeyCode1(KbName(qu1keys{1}))==0 && ...
                    KeyCode1(KbName(qu1keys{2}))==0) && (GetSecs-onset_qu1)<=par.respwi1
                [~,rt1,KeyCode1]=KbCheck;
                 WaitSecs(0.001);
            end
            % computing answer
            if ~isempty(find(KeyCode1,1)) && length(find(KeyCode1))==1
                kbi1=find(KeyCode1);
                switch find(KeyCode1)
                    case KbName(qu1keys{1}), qu1=1;
                    case KbName(qu1keys{2}), qu1=0;
                    otherwise, qu1=NaN;
                end
            else
                kbi1=NaN; 
                qu1=NaN;
            end
            KbReleaseWait;
            
            % Wait between questions to allow PS to stabilize
            Screen(wip(1),'DrawTexture',wip(5)); Screen('Flip',wip(1));
            WaitSecs(par.intrespwi);

            % second question  
            Screen(wip(1),'DrawTexture',wip(7)); [~, onset_qu2]=Screen('Flip',wip(1)); 
            while (KeyCode2(KbName('1!'))==0 && KeyCode2(KbName('2@'))==0 ...
                    && KeyCode2(KbName('3#'))==0 && KeyCode2(KbName('4$'))==0) ...
                    && (GetSecs-onset_qu2)<=par.respwi2
                [~,rt2,KeyCode2]=KbCheck;
                WaitSecs(0.001);
            end
            if ~isempty(find(KeyCode2,1)) && length(find(KeyCode2))==1 
                kbi2=find(KeyCode2);
                switch find(KeyCode2)
                    case KbName('1!'), qu2=1;
                    case KbName('2@'), qu2=2;
                    case KbName('3#'), qu2=3;
                    case KbName('4$'), qu2=4;
                    otherwise, qu2=NaN;
                end
            else
                kbi2=NaN; 
                qu2=NaN;
            end
            KbReleaseWait; 
            
            % Wait after SOJ to allow PS to stabilize
            Screen(wip(1),'DrawTexture',wip(5)); Screen('Flip',wip(1));
            WaitSecs(par.postresp2);

            %%% Deprecated code for response collection at the end of script

            % logging into data arrays
            k=k+1;
            his.kbi(:,k)=[kbi1;kbi2];
            his.qu1(k)=qu1;
            his.qu2(k)=qu2;
            his.ftr(:,k)=[l_rps(i,j);a_rps(i,j)];
            his.timelog(:,k)=[tt-t0;t0-t0;onset_fp-t0;onset_st-t0;offset_st-t0;...
                onset_qu1-t0;rt1-t0;onset_qu2-t0;rt2-t0];
            
            % reinitializing recycled data structures
            KeyCode1=zeros(1,256);
            KeyCode2=zeros(1,256);
            
        end
    end

    Screen('CloseAll');
    Screen('Preference','VisualDebugLevel',oldvdl);
    Screen('LoadNormalizedGammaTable',0,oldgt);
    RestrictKeysForKbCheck([]);
    ShowCursor;
    Priority(0);

    % Save data
    par.timesegment(2,:)=datevec(now);
    if isobject(spo) && strcmp(spo.Status,'open')
        spstr='t_end';
        fprintf(spo,'ET_REM %s\n',spstr); 
    end
    save(par.filename,'par','his');
    
    
    if show_results
        % Converting into score matrix for plotting
        % par.nt x par.nq+2 results array
        resm=zeros(par.nt,par.nq+2);
        resm(1:par.nt,1)=his.ftr(1,1:par.nt);%luminance
        resm(1:par.nt,2)=his.ftr(2,1:par.nt);%additions
        resm(1:par.nt,3)=his.qu1(1:par.nt);
        resm(1:par.nt,4)=his.qu2(1:par.nt);

        scom=zeros(2*par.nq,max(par.nl,par.na));
        for i=1:nt   
            % detection
            scom(1,resm(i,1))=scom(1,resm(i,1))+resm(i,3);  %luminance
            scom(2,resm(i,2))=scom(2,resm(i,2))+resm(i,3);  %additional features
            % second question
            scom(3,resm(i,1))=scom(3,resm(i,1))+resm(i,4);
            scom(4,resm(i,2))=scom(4,resm(i,2))+resm(i,4);
        end

        % Plotting results
        dlm100=scom(1,:)*(100/(par.nt/par.nl));           % respect to luminance values
        slm100=scom(3,:)*(100/(4*par.nt/par.nl));
        dam100=scom(2,:)*(100/(par.nt/par.na));            % respect to additional features 
        sam100=scom(4,:)*(100/(4*par.nt/par.na));

        figure
        subplot(2,2,1)
        plot(lumvec,dlm100)
        axis([0 lumvec(par.nl)+1 0 110])  % axis equal
        set(gca,'TickDir','out')
        xlabel('Target peak luminance in 8bit grayscale space') 
        ylabel(sprintf('Probability of detection ( n=%-d )',par.nt/par.nl))
        title('Detection behavioral psychometric curve')
        subplot(2,2,3)
        plot(lumvec,slm100)
        axis([0 lumvec(par.nl)+1 0 110])  
        set(gca,'TickDir','out')
        xlabel('Target peak luminance in 8bit grayscale space') 
        ylabel(['Across-luminance-average ' qu2mod 's of detection. Scale:[1 2 3 4]'])
        title([qu2mod ' psychometric curve'])
        subplot(2,2,2)
        plot(1:par.na,dam100(1:par.na))
        axis([0 par.na+1 0 110])  
        set(gca,'TickDir','out')
        xlabel('Additional features') 
        ylabel(sprintf('Probability of detection ( n=%-d )',par.nt/par.na))
        title('Detection behavioral psychometric curve')
        subplot(2,2,4)
        plot(1:par.na,sam100(1:par.na))
        axis([0 par.na+1 0 110])  
        set(gca,'TickDir','out')
        xlabel('Additional features') 
        ylabel(['Across-a-average ' qu2mod 's of detection. Scale:[1 2 3 4]'])
        title([qu2mod ' psychometric curve'])
    end
   
catch ME
    
    Screen('CloseAll');
    Screen('Preference','VisualDebugLevel',oldvdl);
    Screen('LoadNormalizedGammaTable',0,oldgt); % repmat((0:1/255:1)',1,3)
    RestrictKeysForKbCheck([]);
    ShowCursor;
    Priority(0);
    
    rethrow(ME);
end



%%% Deprecated code for response  collection
%{
        % presenting stimulus (padded prelaterally with random intervals)
        Screen(wip(1),'DrawTexture',wip(3)); [~, onset_fp]=Screen('Flip',wip(1));
        WaitSecs(par.presi0+par.presi1*rand);
        Screen(wip(1),'DrawTexture',wip(4)); [~, onset_st]=Screen('Flip',wip(1));   
        WaitSecs(par.sdi);
        Screen(wip(1),'DrawTexture',wip(5)); [~, offset_st]=Screen('Flip',wip(1));
        WaitSecs(par.posti);
        
        % detection question
        Screen('CopyWindow',wip(6),wip(1)); [~, onset_detq]=Screen('Flip',wip(1));
        [rt1,KeyCode1]=KbWait([],2,onset_detq+par.respwi1);
        %while find(KeyCode1)~=KbName('Y') && find(KeyCode1)~=KbName('N')...
        %        && GetSecs<onset_detq+par.respwi1
        %    [rt1,KeyCode1]=KbWait([],2,onset_detq+par.respwi1);
        %end
        if ~isempty(find(KeyCode1,1)) && length(find(KeyCode1))==1
            kbi1=find(KeyCode1);
            switch find(KeyCode1)
                case KbName('Y'), whichdet=1;
                case KbName('N'), whichdet=0;
                otherwise, whichdet=NaN;
            end
        else
            kbi1=NaN; 
            whichdet=NaN;
        end
    
        % CR question  
        Screen('CopyWindow',wip(7),wip(1)); [~, onset_crq]=Screen('Flip',wip(1)); % detection CR
        [rt2,KeyCode2]=KbWait([],2,onset_crq+par.respwi2);
        %while find(KeyCode2)~=KbName('1') && find(KeyCode2)~=KbName('2')...
        %        && find(KeyCode2)~=KbName('3') && find(KeyCode2)~=KbName('3')...
        %        && GetSecs<onset_detcrq+par.respwi2
        %   [rt2,KeyCode2]=KbWait([],2,onset_detcrq+par.respwi2);
        %end
        if ~isempty(find(KeyCode2,1)) && length(find(KeyCode2))==1 
            kbi2=find(KeyCode2);
            switch find(KeyCode2)
                case KbName('1!'), qu2=1;
                case KbName('2@'), qu2=2;
                case KbName('3#'), qu2=3;
                case KbName('4$'), qu2=4;
                otherwise, qu2=NaN;
            end
        else
            kbi2=NaN; 
            qu2=NaN;
        end
%}
