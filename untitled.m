% Using edfread
[data, header] = edfread('chb01_01.edf');

info = edfinfo("chb01_01.edf");
disp(info.SignalLabels)