function E3_4discr1(qu2mod,viswon,noam,logdisth,q,pvluam,adpvam,adbglu,spo,lang)
% E3_4: 2AFC orientation(ccw,cw) discrimination task with CR or PAS
% Tropically focused attention: 7(4) values 
%
% nq(2AFC,CR)*ns     *nr= nq*nt
% 02         *(1*3*8)*15= 02*360  uint8 luminance
% Randomize for nt=3*8*15=360=8*9*5=24*5*3 
% 3 blocks of ns*nr=24*5=120


KbName('UnifyKeyNames');
RestrictKeysForKbCheck([KbName('LeftArrow') KbName('RightArrow')...
    KbName('space') KbName('c')...
    KbName('1!') KbName('2@') KbName('3#') KbName('4$')]);
oldvdl=Screen('Preference','VisualDebugLevel',1);
show_results=0;
LoadAcerAF715GammaTable=1;
debugmode=0; 
pilotmode=0;

if viswon
    PsychVideoSwitcher('SwitchMode',0,viswon,0); %high precision luminance mode
end
if pilotmode % Pilot:2 shuffled sequences of 72 trials each: invoke uintrandseq 
    par.nr=6; %nr repetitions of each stimulus
    par.ndisb=2; 
    par.ndist=72;
    load E3_4discr1_rps_pilot2b
else % 3 shuffled sequences of 120 trials each: invoke uintrandseq 
    par.nr=15; %nr repetitions of each stimulus
    par.ndisb=3; 
    par.ndist=120;
    load E3_4discr1_rps
end
par.nrb=par.nr/par.ndisb;  %  nr per block
par.o_rps=o_rps;par.a_rps=a_rps;


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
par.pvluam=pvluam; par.pvluamcl=par.pvluam*par.mgv;
par.lambda=120;            % gabor sine grating spatial frequency 
par.adpvam=adpvam ; par.adpvamcl=par.adpvam*par.mgv; %96 adding extra target contrast set of trials 
par.adbglu=adbglu ; par.adbglucl=par.adbglu*par.mgv;%192 adding extra target average luminance set of trials
par.bglu=0.4; par.bglucl=par.bglu*par.mgv; % stimulus control background luminance
par.noam=noam; par.noamcl=par.noam*par.mgv;      % noise amplitude
gamma=0.5;grain=0.01;
par.disth=10*10.^logdisth;
ivq=find(QuestP(q,-1:grain:1)>=gamma+0.05,1,'first');
vq=QuestP(q,ivq*grain-1);  % virgintile=20-quantile
logazvec=logdisth*ones(1,4)+vq*[-inf,-3/5,1/5,1];
azvec=10*10.^logazvec;  % target azimuth values
%The mean and variance are always specified as if the image were of class double in the range [0 1]. 
%If the input image is of class uint8 or uint16, the imnoise function converts the image to double, 
%adds noise according to the specified type, and then converts the noisy image back to the same class as the input
% Defining paremeters: trial sequence
par.orivec=90*ones(1,8)+[-1*azvec(4),-1*azvec(3),-1*azvec(2),0,0,azvec(2),azvec(3),azvec(4)]; % target orientation values
par.no=length(par.orivec);          % number of target peak luminance values  
par.na=3;                        % control, additional average luminance series, and additional spatial frequency series  
par.ns=1*par.na*par.no;                    % different stimuli types
par.nt=par.ndist*par.ndisb;                    % total number of trials nt=nr*ns=par.ndist*par.ndisb
par.nq=2;                         % number of questions per trial   
par.presti=1;           % prestimulus interval: baseline for PS
par.presi0=.3; par.presi1=.3;         % prestimulus interval: randomizing onset time
par.sdi=.3;                   % stimulus display time in ms
par.posti=0;                     % poststimulus interval: allow immediate foj
par.intrespwi=1;             % gap to measure PS before soj
par.respwi1=inf;                   % response window for detection
par.respwi2=inf;                   % response window for CR
par.postresp2=1;                  % window for measuring PS after soj
%ITI=presi0+presi1*rand+sdi+posti+respwi=300+300+300:300+300+300+300+respwi=900:2400+respwi
%The small-angle approximation holds for angles up to ~10degrees  
%The angles at which the relative error exceeds 1% is sin ??????ｽ?????ｽ??????ｽ?????ｽ???????ｽ?????ｽ??????ｽ?????ｽat about 0.244rad=13.98deg


% Stacking stimuli frames in buffer cells
[~,par.l,c,par.csrad,par.sigma]=E3gabor_noise(par.bglu,par.noam,par.lambda,45,5);       % extracting image of size l
EM64=par.bglu*ones(par.l);         % empty matrix 
FM64=EM64;
FM64(c-3:c+3,c)=par.fplu;FM64(c,c-3:c+3)=par.fplu;     % crosshead FP 
BC64=cell(par.no,par.na);             % cell holding nc images in buffer  
vadbglu=[1 adbglu/par.bglu 1];        % adding extra target average luminance set of trials
vadpvam=[1 1 adpvam/pvluam];      % adding extra target contrast set of trials
for i=1:par.no 
    for j=1:par.na
        BC64{i,j}=E3gabor_noise(par.bglu*vadbglu(j),par.bglu,par.noam,...
            par.lambda,par.orivec(i),par.pvluam*vadpvam(j)); % combining OR with additives  
    end
end

if viswon
    EM14=PsychVideoSwitcher('MapLuminanceToRGB',EM64,par.btrr);
    FM14=PsychVideoSwitcher('MapLuminanceToRGB',FM64,par.btrr);
    BC14=BC64;
    for i=1:par.no 
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
    for i=1:par.no 
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
if strcmp(lang,'EN')
    instr{1}='Orientation discrimination task';
    instr{2}=' '; %This task will measure your perceptual sensitivity, and your subjective experience.';
    instr{3}='Gratings with different orientations will be presented at the center of the screen.';
    instr{4}='By selecting the corresponding keys, you will have to report:';
    instr{5}='Orientation discrimination (ccw,cw)';
    instr{6}='[Left]: Counterclockwise-tilted   [Right]: Clockwise-tilted';
    if strcmp(qu2mod,'CR')
        instr{7}='Confidence in your previous decision';
        instr{8}='[1]: not confident at all  [2]: a little confident  [3]: quite confident  [4]: absolutely confident';
    elseif strcmp(qu2mod,'PAS')
        instr{7}='Visibility of the stimulus orientation';
        instr{8}='[1]: seen not tilted  [2]: seen vaguely tilted  [3]: seen almost clearly tilted  [4]: seen absolutely clearly tilted';
    else 
        error('KuvikException:InvalidArgument','Unrecognized function argument.');
    end
    instr{9}='IMPORTANT: try to use the {1 2 3 4} scale as exhaustively as possible';
    instr{10}='IMPORTANT: answer as fast as possible the first question after the stimulus, but give priority to response accuracy';
    instr{11}='You will receive your performance score at the end of the experiment';
    instr{12}='Press spacebar when ready';
    instr{13}='{<-  ->}';
    instr{14}='{1 2 3 4}';
    switch OSName
        case 'Linux', fns='-misc-fixed-medium-r-normal--18-120-100-100-c-90-iso8859-15';
        case 'Windows', fns='Helvetica';
        otherwise, fns='Helvetica';
    end
elseif strcmp(lang,'JP')
    instr{1}='傾き識別課題';
    instr{2}=' ';
    instr{3}='異なる傾きの縞模様が画面の中央に表示されます。';
    instr{4}='対応するキーを押すことで、次の模様の特徴を報告してください。';
    instr{5}='傾き識別';
    instr{6}='[<-]:反時計回り　　[->]:時計回り';
    if strcmp(qu2mod,'CR')
        instr{7}='直前の回答に対するあなたの確信度';
        instr{8}='[1]:全く自信がない  [2]:少し自信がある  [3]:かなり自信がある  [4]:完全に自信がある';
    elseif strcmp(qu2mod,'PAS')
        instr{7}='縞模様の見え方';
        instr{8}='[1]:傾いて見えなかった  [2]:わずかに傾いて見えた  [3]:少し傾いて見えた  [4]:完全に傾いて見えた';
    else 
        error('KuvikException:InvalidArgument','Unrecognized function argument.');
    end
    instr{9}='重要：確信度を回答する際、課題に合わせてなるべく{1 2 3 4}スケールのすべての要素を使用して下さい。';
    instr{10}='重要：一番目の質問にできるだけ早く答えて下さい。ただし、早さより正確さを優先してください。';
    instr{11}='実験後、課題のスコアをお知らせします。';
    instr{12}='準備ができたら、スペースキーを押して始めてください。';
    instr{13}='{<-  ->}';
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
par.filename=['E3_4discr1' qu2mod '_' dateISO8601];
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
Screen(wip(8),'TextSize',20); bourec=Screen('TextBounds', wip(8),instr{2});
Screen(wip(8),'DrawText',instr{2},(par.scrsiz(3)-bourec(3))/2,.1*par.scrsiz(4));
bourec= Screen('TextBounds', wip(8),instr{3});
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
bourec= Screen('TextBounds', wip(8),instr{9});
Screen(wip(8),'DrawText',instr{9},(par.scrsiz(3)-bourec(3))/2,.65*par.scrsiz(4));
bourec= Screen('TextBounds', wip(8),instr{10});
Screen(wip(8),'DrawText',instr{10},(par.scrsiz(3)-bourec(3))/2,.7*par.scrsiz(4));
bourec= Screen('TextBounds', wip(8),instr{11});
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


%%%%%%%%%%%%%%%%%
% task
%%%%%%%%%%%%%%%%%

HideCursor;
par.prilev=MaxPriority(wip(1));
Priority(par.prilev);

Screen(wip(1),'DrawTexture',wip(8)); 
[~, tt]=Screen('Flip',wip(1));

[t0,KeyCode]=KbWait([],2);
while ~KeyCode(KbName('space'))==1, [t0,KeyCode]=KbWait([],2); end
if isobject(spo) && strcmp(spo.Status,'open')
    spstr='t0';
    fprintf(spo,'ET_REM %s\n',spstr); 
end
    
k=0;

for j=1:par.ndisb
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
       
    for i=1:par.ndist
        % Putting stimulus in buffer 4
        wip(4)=Screen(wip(1),'MakeTexture',BC{o_rps(i,j),a_rps(i,j)});
                          
        % Baseline window
        WaitSecs(par.presti);
            
        % presenting stimulus (padded prelaterally with random intervals)
        Screen(wip(1),'DrawTexture',wip(3)); [~, onset_fp]=Screen('Flip',wip(1));
        WaitSecs(par.presi0+par.presi1*rand);
        Screen(wip(1),'DrawTexture',wip(4)); [~, onset_st]=Screen('Flip',wip(1));   
        if isobject(spo) && strcmp(spo.Status,'open') 
            spstr='onset_st';%['E3_4discr1',qu2mod,' onset_st-',num2str(o_rps(i,j)),'-',num2str(a_rps(i,j))];
            fprintf(spo,'ET_REM %s\n',spstr);  
        end
        % show stimulus until until sdi elapses
        while (GetSecs-onset_st)<=par.sdi
            % poll for a response: subjects can respond before stimulus terminates
            if (KeyCode1(KbName('LeftArrow'))==1 || KeyCode1(KbName('RightArrow'))==1)
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
                
        % orientation discrimination question: loop until valid key is pressed or reswi terminates
        Screen(wip(1),'DrawTexture',wip(6)); [~, onset_qu1]=Screen('Flip',wip(1));
        while (KeyCode1(KbName('LeftArrow'))==0 && KeyCode1(KbName('RightArrow'))==0) ...
                && (GetSecs-onset_qu1)<=par.respwi1
             [~,rt1,KeyCode1]=KbCheck;
             WaitSecs(0.001);
        end
        % computing accuracy
        if ~isempty(find(KeyCode1,1)) && length(find(KeyCode1))==1
           kbi1=find(KeyCode1);
           switch find(KeyCode1)
               case {KbName('LeftArrow'),KbName('RightArrow')}, 
                   qu1=and(kbi1==KbName('LeftArrow'),o_rps(i,j)>4.5)||...
                       and(kbi1==KbName('RightArrow'),o_rps(i,j)<4.5);     %------WARNING: >4.5 is not exactly accuracy, this is an approximation ------%
               otherwise, qu1=nan;
           end
        else
            kbi1=nan;
            qu1=nan;
        end
        KbReleaseWait; 
        
        % Wait between questions to allow PS to stabilize
        Screen(wip(1),'DrawTexture',wip(5)); Screen('Flip',wip(1));
        WaitSecs(par.intrespwi);

        % second question  
        Screen(wip(1),'DrawTexture',wip(7)); [~, onset_qu2]=Screen('Flip',wip(1));
        while ( KeyCode2(KbName('1!'))==0 && KeyCode2(KbName('2@'))==0 ...
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
                otherwise, qu2=nan;
            end
        else
            kbi2=NaN; 
            qu2=NaN;
        end
        KbReleaseWait; 
        
        % Wait after SOJ to allow PS to stabilize
        Screen(wip(1),'DrawTexture',wip(5)); Screen('Flip',wip(1));
        WaitSecs(par.postresp2);
            
        % logging into data arrays
        k=k+1;
        his.kbi(:,k)=[kbi1;kbi2];
        his.qu1(k)=qu1;
        his.qu2(k)=qu2;
        his.ftr(:,k)=[o_rps(i,j);a_rps(i,j)];
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
    resm(1:par.nt,1)=his.ftr(1,1:par.nt); %orientation
    resm(1:par.nt,2)=his.ftr(2,1:par.nt); %additions
    resm(1:par.nt,3)=his.qu1(1:par.nt);
    resm(1:par.nt,4)=his.qu2(1:par.nt);

    scom=zeros(2*nq,max(par.no,par.na));
    for i=1:par.nt   
        % 2AFC discrimination
        scom(1,resm(i,1))=scom(1,resm(i,1))+resm(i,3);  %orientation
        scom(2,resm(i,2))=scom(2,resm(i,2))+resm(i,3);  %additional features
        % 2AFC discrimination CR
        scom(3,resm(i,1))=scom(3,resm(i,1))+resm(i,4);
        scom(4,resm(i,2))=scom(4,resm(i,2))+resm(i,4);
    end

    %
    % Plotting results
    dom100=scom(1,:)*(100/(par.nt/par.no));           % respect to orientation
    som100=scom(3,:)*(100/(4*par.nt/par.no));
    dam100=scom(2,:)*(100/(par.nt/par.na));      % respect to additional features 
    sam100=scom(4,:)*(100/(4*par.nt/par.na));

    figure
    subplot(2,2,1)
    plot(orivec,dom100)
    axis([0 lumvec(par.no)+1 0 110])  % axis equal
    set(gca,'TickDir','out')
    xlabel('Stimulus orientation in degrees') 
    ylabel(sprintf('Probability of discrimination ( n=%-d )',par.nt/par.no))
    title('2AFC discrimination behavioral psychometric curve')
    subplot(2,2,3)
    plot(orivec,som100)
    axis([0 lumvec(par.no)+1 0 110])  
    set(gca,'TickDir','out')
    xlabel('Stimulus orientation in degrees') 
    ylabel(['Across-orientation-average ' qu2mod 's of orientation discrimination. Scale:[1 2 3 4]'])
    title([qu2mod ' psychometric curve'])
    subplot(2,2,2)
    plot(1:par.na,dam100(1:par.na))
    axis([0 par.na+1 0 110])  
    set(gca,'TickDir','out')
    xlabel('Additional features') 
    ylabel(sprintf('Probability of discrimination ( n=%-d )',par.nt/par.na))
    title('2AFC discrimination behavioral psychometric curve')
    subplot(2,2,4)
    plot(1:par.na,sam100(1:par.na))
    axis([0 par.na+1 0 110])  
    set(gca,'TickDir','out')
    xlabel('Additional features') 
    ylabel(['Across-a-average ' qu2mod 's of detection. Scale:[1 2 3 4]'])
    title([qu2mod ' psychometric curve'])
end





