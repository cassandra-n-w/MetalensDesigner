
% note can use uigetdir to get a path like this. quite handy
folderpath = 'C:\Users\Cassie\OneDrive\Drive\Documents\My Work\Research\metamaterial_lens\Measurements\480GHz Measurements\Cassie Measurements 3-29-23';

files = dir(fullfile(folderpath, "*.dat"));


for j = 1:length(files)

    name = files(j).name;

    Measures(j) = Meas(folderpath, name);
end

Ycut = 4.05;
Xcut = 11.30;

freq = 500e9;

Ycut_Xs = [];
Xcut_Ys = [];

Ycut_Db = [];
Xcut_Db = [];

for i = 1:length(Measures)
    meas = Measures(i);

    if meas.X == Xcut
        Xcut_Ys = [Xcut_Ys, meas.Y];
        Xcut_Db = [ Xcut_Db, meas.select_data(freq)];
    end

    if meas.Y == Ycut
        Ycut_Xs = [Ycut_Xs, meas.X];
        Ycut_Db = [Ycut_Db, meas.select_data(freq)];
    end
end

hold off;
Ycut_Xs = Ycut_Xs - 11.3;
Xcut_Ys = Xcut_Ys - 4.5;
plot(Ycut_Xs, Ycut_Db, ".", Xcut_Ys, Xcut_Db, ".")

xlabel("Horn Position (mm)");
ylabel("S21 (dB)");
ylim([-60, 0])
title(strcat("2f - 2f Horn Coupling Lens Test, ", num2str(freq/1e9), "GHz"))

hold on;
if exist("xcut", "var")
    idx = interp1(freqs, 1:length(freqs),freq,'nearest'); 
    meas_dB = mag2db(abs(xcut(:,:,idx)));
    horn_loss = 1 - 0.0947;
    horn_loss_dB = 10*log10(horn_loss);
    meas_dB = meas_dB + 2*horn_loss_dB;
    plot(xscan, meas_dB, yscan, meas_dB - 1.9);
    legend(["X cut", "Y cut", "Simulated", "Simulated - 1.9dB"]);
else
    legend(["X cut", "Y cut"]);
end
hold off;




