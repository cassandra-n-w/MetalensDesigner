
% note can use uigetdir to get a path like this. quite handy
folderpath = 'C:\Users\Cassie\OneDrive\Drive\Documents\My Work\Research\metamaterial_lens\Measurements\480GHz Measurements\Cassie Measurements 3-29-23';

files = dir(fullfile(folderpath, "*.dat"));


for j = 1:length(files)

    name = files(j).name;

    Measures(j) = Meas(folderpath, name);
end

Ycut = 4.05;
Xcut = 11.30;

freq = 480e9;

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
title("2f - 2f Horn Coupling Lens Test, 480GHz")

hold on;
if exist("xcut")
    plot(xscan, mag2db(abs(xcut(:,:,4))), yscan, mag2db(abs(ycut(:,:,4))) - 2.8);
    legend(["X cut", "Y cut", "Simulated", "Simulated - 2.8dB"]);
else
    legend(["X cut", "Y cut"]);
end




