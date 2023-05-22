
% note can use uigetdir to get a path like this. quite handy
folderpath = 'C:\Users\Cassie\OneDrive\Drive\Documents\My Work\Research\metamaterial_lens\Measurements\480GHz Measurements\semiauto_measurement_1_5-17-23';

files = dir(fullfile(folderpath, "*.mat"));

Xspacing = 0.2; % 0.2mm spacing was chosen between x values in auto scan
for j = 1:length(files)

    name = files(j).name;

    Measures(j) = SemiautoMeas(folderpath, name, Xspacing);
end

[~, sort_ind] = sort([Measures.Y]);
Measures = Measures(sort_ind);

Ycut = 4.05;
Xcut = 11.30;

freq = 480e9;

Ycut_Xs = [];
Xcut_Ys = [];

Ycut_Db = [];
Xcut_Db = [];

fulldata = [];
xrange = [];
yrange = [];

for i = 1:length(Measures)
    meas = Measures(i);

    fulldata = [fulldata meas.select_data(freq)];
    xrange = meas.X;
    yrange = [yrange meas.Y];
end


%imagesc(xrange, yrange, flipud(transpose(abs(fulldata))));
%title("2D lens coupling scan, linear |S21| measurement");
imagesc(xrange, yrange, flipud(transpose(mag2db(abs(fulldata)))));
title("2D lens coupling scan, power measurement in dB");
clim([-50,0])

xlabel("X position (mm)");
ylabel("Y position (mm)");

axis equal
xlim([min(xrange), max(xrange)]);

colorbar



