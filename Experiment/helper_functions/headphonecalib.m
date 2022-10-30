function [calibrationleft, calibrationright] = headphonecalib(freq)
% give matlab sound level meter data
frequency = [500 750 1000 1500 2000 3000 4000 5000];
leftamplitude = [0.5 0.4 0.35 0.35 0.4 0.6 0.7 0.2];
rightamplitude = [0.5 0.45 0.4 0.35 0.4 0.7 1.1 0.2];

%interpolate all frequencies using recorded data 
xx = 1:1:5000;
interpolatedleft = spline(frequency,leftamplitude,xx);
interpolatedright = spline(frequency, rightamplitude, xx);

%find amplitude for frequency 
calibrationleft = interpolatedleft(freq);
calibrationright = interpolatedright(freq);
end 