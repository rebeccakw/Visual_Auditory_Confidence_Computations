if strcmp(P.stim_type, 'audio')
    beep.pahandle = PsychPortAudio('Open', [], 1, 1, P.samplingrate, P.nrchannels);
end

for category = 1 : 2
    DrawFormattedText(scr.win, sprintf('Examples of Category %i stimuli:', category), 'center', scr.cy-60, color.wt);
    Screen('Flip', scr.win);
    WaitSecs(2.5);
    
    for i = 1:nDemoTrials
        
        Screen('DrawTexture', scr.win, scr.cross);
        Screen('Flip', scr.win);
        WaitSecs(Demo.t.betwtrials/1000);
        
       
        if strcmp(P.stim_type, 'audio')
          stim.ort = round(stimulus_orientations(Test.category_params, category, 1, category_type),0);
        else
          stim.ort = stimulus_orientations(Test.category_params, category, 1, category_type);   
        end
       
        if strcmp(P.stim_type, 'grate')
            stim.phase = 360*rand;
            grate(P, scr, Demo.t, stim);
        else
            audio(P, scr, stim, beep); 
            WaitSecs(0.75);
        end
    end
    
    flip_key_flip(scr,'continue','center',color,new_subject);
    
end

if strcmp(P.stim_type, 'audio')
PsychPortAudio('Close', beep.pahandle)
end 