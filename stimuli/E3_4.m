function E3_4(qu2mod,viswon,noam,spo,lang)
% E3_4:  Accuracy and awareness as sigmoid functions of stimulus feature perceptibility
% Measures: PMF, pupil size
% Probing: interactions in 
%    Type2PMF x(intra) lowlevelfeatures x(inter)Type2Mode | processingdepth
% Methods: 4-point scales will be used for both CR and PAS
% #trials: nt=720+360+270=1350

% E3_4pilot:
% #trials: nt=120+288+144+108=120+540

if exist('rng','builtin'), rng('default'); else rand('state',0); end

% Parameters
grain=0.01; gamma=0.5;
bglu=0.4;
lambda=120;
%noam=0;
% Condition parameter: qu2mod


%% E3.4.0: Y/N detection task
%viswon=0;noam=0;bglu=0.4; qu2mod='CR';

pThreshold0=0.75;
nt0=50;

[logdethr0,q0]=E3calib_det(pThreshold0,viswon,noam,nt0,spo,lang); % >logdethr
EEndScr(0,lang)

% FOR FIGURE MAKING
if 1, logdethr=log10(0.001/0.4), q0=zeros(1,100),end

% Constant stimuli, randomized blocks
% Monotropical(vertical) throughout all the trials: uniformize attentional effects
% One relevant feature: contrast
% Two additional irrelevant features: spatial frequency and background luminance
adbglu0=11/10*bglu;     % background luminance
adsfr0=4/5*lambda;  % spatial frequency   
% used luminance range maximum: 0.4*1.05+0.2+dethr*10=0.62+dethr*10
% used luminance range minimum: 0.4-0.2-dethr*10=0.2-dethr*10
ivq0=find(QuestP(q0,-1:grain:1)>=gamma+0.05,1,'first');
vq0=QuestP(q0,ivq0*grain-1);
maxpvlu0=bglu*10^(logdethr0+vq0*21/20); 
if ~isempty(logdethr0) && maxpvlu0>0.2       %logdethr+vq*21/20 > -.301  
    error('KuvikException:ExceededRange',...
        'Out of luminance range. Noise distribution will be rectified in the main task.');
end
% Issue exceeded luminance range error if logdethr+vq*21/20 > -.301     

% Y/N detection task with CR or PAS
% nq(Y/N,CR)*ns     *nr= nq*nt
% 02        *(1*3*8)*30= 02*720n
% randomized blocking: 6 blocks of 120
E3_4detec(qu2mod,viswon,noam,logdethr0,q0,adbglu0,adsfr0,spo,lang);               
EEndScr(0,lang)



%% E3.4.1: 2AFC orientation(ccw,cw) discrimination task 1 
% Constant stimuli, randomized blocks
% Tropically focused attention: 7(4) values 
%viswon=0;noam=0.0;bglu=0.4; qu2mod='CR';
Threshold11=0.75; 
nt11=35;

% Not necessary: better to fix to logdethr1=-1
%Threshold10=0.75;
%logdethr1=E3calib_det(Threshold10,viswon,noam,nt10); % >logdethr
%EEndScr(0,lang)
logdethr1=-1;

% One relevant feature: orientation
% Two additional irrelevant features: background luminance and target contrast
pvluam1=bglu*10^logdethr1;      %luminance of the gabor peak 
adpvam1=7/6*bglu*10^logdethr1;      %added luminance of gabor peak(contrast)
adbglu1=11/10*bglu;                        %added background luminance
% used luminance range maximum: 0.4*1.1+0.2+dethr*2=0.62+dethr*2
% used luminance range minimum: 0.4-0.2-dethr*2=0.2-dethr*2
maxpvlu1=2*bglu*10^logdethr1; 
if ~isempty(logdethr1) && maxpvlu1>0.2       %logdethr > -.6021   
    error('KuvikException:ExceededRange',...
        'Out of luminance range. Noise distribution will be rectified in the main task.');
end
% Issue exceeded luminance range warning if logdethr > -.6021

[logdisth1,~,q1]=E3calib_dis1(Threshold11,viswon,noam,pvluam1,nt11,spo,lang); % <pvluam 
EEndScr(0,lang)

% 2AFC orientation(left,right) discrimination task with CR or PAS
% nq(2AFC,CR)*ns     *nr= nq*nt
% 02         *(1*3*8)*15= 02*360
% randomized blocking: 3 blocks of 120
E3_4discr1(qu2mod,viswon,noam,logdisth1,q1,pvluam1,adpvam1,adbglu1,spo,lang);
%save(['E3_4discr1CR' cglog_fn(7:26)],'scom')%saving results to mat file
EEndScr(0,lang)



%% E3.4.2: 4AFC orientation(0,45,90,135) discrimination task 2  

% Constant stimuli, randomized blocks
% Tropically scattered attention: 4 values 
% IMPORTANT: tell the subject to press always Y
pThreshold2=0.625;
nt2=35;

[logdisth2,~,q2]=E3calib_dis2(pThreshold2,viswon,noam,nt2,spo,lang);
EEndScr(0,lang)

% One relevant feature: contrast
% Two additional irrelevant features: spatial frequency and background luminance
adbglu2=11/10*bglu; 
adsfr2=4/5*lambda;    
% used luminance range maximum: 0.4*1.1+0.2+dis2th*3=0.62+dis2th*3
% used luminance range minimum: 0.4-0.2-dis2th*3=0.2-dis2th*3
ivq2=find(QuestP(q2,-1:grain:1)>=gamma+0.05,1,'first');
vq2=QuestP(q2,ivq2*grain-1);
maxpvlu2=bglu*10^(logdisth2+vq2*1); 
if ~isempty(logdisth2) && maxpvlu2>0.2       %logdisth+vq*1 > -.301   
    error('KuvikException:ExceededRange',...
        'Out of luminance range. Noise distribution will be rectified in the main task.');
end
% Issue exceeded luminance range error if logdisth+vq*3/5 > -.301   

% 4AFC orientation(ccw,cw) discrimination task with CR or PAS
% nq(4AFC,CR)*ns     *nr= nq*nt
% 02         *(1*3*3)*30= 02*270
% randomized blocking: 3 blocks of 90
E3_4discr2(qu2mod,viswon,noam,logdisth2,q2,adbglu2,adsfr2,spo,lang);
%save(['E3_4discr2CR' cglog_fn(7:26)],'scom')%saving results to mat file
EEndScr(0,lang)
