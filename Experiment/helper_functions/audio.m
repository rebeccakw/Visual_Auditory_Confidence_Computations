function [truebeeplength] = audio(P, scr, stim, beep)
%flip to grey screen 
t0 = Screen('Flip', scr.win);

% get values to equalise headphones
[left, right] = headphonecalib(stim.ort);

% Make a beep which we will play back to the user
leftBeep = left*MakeBeep(stim.ort, P.beepLengthSecs, P.samplingrate);
rightBeep = right*MakeBeep(stim.ort, P.beepLengthSecs, P.samplingrate);

%Apply ramp 
beeponly = P.beepLengthSecs-2*P.rampLengthSecs;
ramp = [linspace(0,1,P.rampLengthSecs*P.samplingrate), repelem(1,(round(beeponly*P.samplingrate,0))), linspace(1,0,P.rampLengthSecs*P.samplingrate)];
leftBeep = ramp.*leftBeep;
rightBeep = ramp.*rightBeep;

%apply intensity correction 
leftBeep = (stim.sigma).*leftBeep; %adjustment to tone not in log10 units
rightBeep = (10^intensity).*rightBeep;


% Fill the audio playback buffer with the audio data, doubled for stereo
% presentation
PsychPortAudio('FillBuffer', beep.pahandle, [rightBeep; leftBeep]);

% Start audio playback
PsychPortAudio('Start', beep.pahandle, 1, 0, 1);
 
% Wait for stop of playback
[actualStartTime, ~, ~, estStopTime] = PsychPortAudio('Stop', beep.pahandle, 1, 1);

truebeeplength = estStopTime - actualStartTime;

end