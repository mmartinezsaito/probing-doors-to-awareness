
p=20;
debugmode=0;
viswon=0;

switch debugmode
    case 0, winsiz=[];
    case 1, winsiz=[+10 +10 par.scrsiz(3)/2-10 par.scrsiz(4)/2-10];
    case 2, PsychDebugWindowConfiguration([],.9), winsiz=[];
end      

wip=Screen(0,'OpenWindow',[],winsiz);

Screen(wip,'FillRect',130);
      
                                                                   
for i=0:p
                                                                                                   
    I=ones(700)*((129+i/p)/255);
    
         
    if viswon
        I=PsychVideoSwitcher('MapLuminanceToRGB',I,126.3);        
    else                       
        I=uint8(I*255);    
    end
    
                                           
    wipt=Screen(wip,'MakeTexture',I);  
    
    
    KbWait([],2);                     
    
    
    Screen(wip,'DrawTexture',wipt);
    Screen('Flip',wip);
    
  
end

sca
