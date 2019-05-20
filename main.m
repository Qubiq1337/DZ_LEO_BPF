clear; 
Fs = 800;
ts = 0:1.0/Fs:1-0.0000001;
freq = 10;
edit input.txt;
x = sin(2*pi*freq*ts) + sin(20*pi*freq*ts) + sin(40*pi*freq*ts);
fileID = fopen('input.txt', 'w');
fprintf(fileID,'%f\n', x);
fprintf(fileID,'%i', Fs);
fclose(fileID);

