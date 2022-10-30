function countdown(scr,color,x,y)
for i=1:scr.countdown_time+1;
    WaitSecs(1);
    Screen('FillRect',scr.win,color.bg,[x-4 y-40 x+1.5*scr.fontsize y+0.5*scr.fontsize]) %timer background
    DrawFormattedText(scr.win,[num2str(scr.countdown_time+1-i) '  '],x,y,color.wt);
    Screen('Flip',scr.win,[],1); % flip to screen without clearing
end
end
