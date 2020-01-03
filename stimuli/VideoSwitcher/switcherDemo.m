% This demo shows how to use video switcher to achieve low contrast gabor.
% Run it and you follow the screen instruction.

% 11/2006   wrote it (xl)
% 03/2011   use PsychImaging (xl)
 
function switcherDemo(whichScreen)
% display parameters
if nargin<1, whichScreen=max(Screen('screens')); end
isBox=1;        % change to 0 if you are using a video switcher card
btrr=126.3;     % blue to red ratio of video switcher
ppd=84;         % pixels per degree based on monitor and viewing distance
gamma=2;        % monitor gamma

% stimulus parameter
sf=1;           % cycle per degree
angle=0 ;       % initial orientation
con=0.2;        % initial contrast
radius=2;       % image radius in degree
back=0.5;       % normalized background
useTrigger=0;   % enable trigger in pipline

% keys to control stimulus
KbName('UnifyKeyNames');
keys={'UpArrow' 'DownArrow' 'LeftArrow' 'RightArrow' 'Escape'};

try % avoid dead screen in case of error
    PsychVideoSwitcher('SwitchMode', whichScreen, 1, isBox); % swtich into gray scale mode
 
    % PsychImaging pipelinwe setup
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32Bit');
    PsychImaging('AddTask', 'General', 'EnableVideoSwitcherSimpleLuminanceOutput', btrr, useTrigger);
    PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
    [w r]= PsychImaging('OpenWindow', whichScreen, back); %#ok
    HideCursor;
    PsychColorCorrection('SetEncodingGamma', w, 1/gamma);
    
    % set larger text for instruction
    Screen('TextFont',w,'Times');
    Screen('TextSize',w,24);
    
    m=round(ppd*radius*2); % image size in pixels
    sig=m/6;               % gabor sigma
    sf=sf/ppd;             % cycles per pixel
    
    % create the only texture for gabor
    tex=CreateProceduralGabor(w,m,m,0,[1 1 1 0]*back,1,0.5);

    while 1  % loop till ESC is pressed
        
        % print instruction and result
        x=Screen('DrawText',w,'Up/Down arrow keys to change contrast by 10% (now ',20,10);
        Screen('DrawText',w,num2str(con,'%.2g).'),x,10);
        x=Screen('DrawText',w,'Left/Right arrow keys to rotate grating by 5 degrees (now ',20,50);
        Screen('DrawText',w,num2str(angle,'%g).'),x,50);
        Screen('DrawText',w,'ESC to quit.',20,90);
        
        % draw gabor into buffer, all parameters are passed on the fly
        Screen('DrawTexture',w,tex,[],[],angle,[],[],[],[],2,[90 sf sig con 1 0 0 0]);
        
        % this if block is for trigger
        if useTrigger
            % set trigger at vertical center, for only next flip
            PsychVideoSwitcher('SetTrigger', w, r(4)/2, 1); %#ok
            
            % show stimulus and send trigger
            Screen('Flip',w, 0, 1); % don't clear buffer, so we can turn off trigger
            
            % Disable trigger at next flip. This is needed here since we
            % won't do frame by frame change. If you flip each frame, then
            % you will SetTrigger for one flip, and don't need above extra
            % non-clear flip and next disable line.
            PsychVideoSwitcher('SetTrigger', w, []);
        end

        Screen('Flip',w); % show stimulus and turn off trigger if used
        
        % detect keys to make change to gabor
        KbReleaseWait; % avoid repeated key detection
        switch WaitTill(keys)
            case 'UpArrow'
                con=con*1.1; con=min(con,1);
            case 'DownArrow'
                con=con/1.1; con=max(con,0.0005);
            case 'LeftArrow'
                angle=angle+5; angle=mod(angle,360);
            case 'RightArrow'
                angle=angle-5; angle=mod(angle,360);
            otherwise, break;    % esc
        end
    end
    
catch me 
end

Screen('CloseAll');	% close screen and show cursor
PsychVideoSwitcher('SwitchMode', whichScreen, 0, isBox); % switch back to color mode

if exist('me','var'), rethrow(me); end % show error message if any

