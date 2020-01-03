function EEndScr(t,lang)

KbName('UnifyKeyNames');
StuckKeysDisabler;
oldvdl=Screen('Preference','VisualDebugLevel',1);
HideCursor;
LoadAcerAF715GammaTable=1;

bglu=0.4;
texlu=0.6;
mgv=2^8-1;

scrsiz=get(0,'ScreenSize');  %[scrsiz(3) scrsiz(4)]=Screen('WindowSize', 0);
wip=Screen(0,'OpenWindow',[],[]);
Screen(wip,'FillRect',bglu*mgv);
oldgt=Screen('ReadNormalizedGammaTable',0);
if LoadAcerAF715GammaTable
    load AcerAF715; graygammatable=gammaTable1;
    Screen('LoadNormalizedGammaTable',0, graygammatable*[1 1 1]);
end
      
wfns='Helvetica';
lfns='-misc-fixed-medium-r-normal--18-120-100-100-c-90-iso8859-15';
switch OSName
    case 'Linux',fns=lfns;
    case 'Windows', fns=wfns;
    otherwise, fns=wfns;
end
if strcmp(lang,'EN')
    if ~t
        instr{1}='End of this task';
        instr{2}='You can take a rest';
        instr{3}='Press any key when ready for the next task';
    elseif t
        instr{1}='End of experiment';
        instr{2}='Thank you!';
        instr{3}='Press any key to exit';
    end
    switch OSName
        case 'Linux', fns='-misc-fixed-medium-r-normal--18-120-100-100-c-90-iso8859-15';
        case 'Windows', fns='Helvetica';
        otherwise, fns='Helvetica';
    end
elseif strcmp(lang,'JP')
        if ~t
        instr{1}='本課題終了。続きがあります。';
        instr{2}='休憩しても構いません';
        instr{3}='次の課題に進んでもよければ任意のキーを押して始めてください';
    elseif t
        instr{1}='実験終了';
        instr{2}='ありがとうございました';
        instr{3}='任意のキーを押してください';
        end
    instr=cellfun(@double,instr,'uniformOutput',false);
    fns='TakaoExGothic';
    if IsLinux, fns='-:lang=ja'; end
end

wipt=Screen(0,'OpenOffscreenWindow',bglu*mgv);
Screen(wipt,'TextColor',texlu*mgv*[1 1 1]);
Screen(wipt,'TextFont',fns);
Screen(wipt,'TextSize',30); bourec=Screen('TextBounds',wipt,instr{1});
Screen(wipt,'DrawText',instr{1},(scrsiz(3)-bourec(3))/2,.4*scrsiz(4));
Screen(wipt,'TextSize',20); bourec=Screen('TextBounds',wipt,instr{2});
Screen(wipt,'DrawText',instr{2},(scrsiz(3)-bourec(3))/2,.6*scrsiz(4));
Screen(wipt,'TextSize',20); bourec=Screen('TextBounds',wipt,instr{3});
Screen(wipt,'DrawText',instr{3},(scrsiz(3)-bourec(3))/2,.9*scrsiz(4));

Screen(wip,'DrawTexture',wipt); 
Screen('Flip',wip);

KbWait([],2);                     

Screen('Close',wip);
Screen('LoadNormalizedGammaTable',0,oldgt);   % repmat((0:1/255:1)',1,3)
Screen('Preference','VisualDebugLevel',oldvdl);
ShowCursor;
