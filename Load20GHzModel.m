load('matlab_results.mat');

f0 = f;

freqs = linspace(f*0.9, f*1.1, 30);

el = Element('C:\\Users\\Cassie\\Google Drive\\Dropbox\\Documents\\My Work\\Research\\metamaterial_lens\\code\\HFSS Models\\RO3003 touchstone');
di = Dielectric(3.00 * ones(size(freqs)), 0.001, freqs, 0.5, 1.5);

proto = ProtoTL(el, di);

lensmodel = Lens(dim(1), layer_thickness*ones(11, 1),  zf, g, proto);
%lensmodel = Lens(dim(1), 0.8*layer_thickness*ones(11, 1),  zf, g, proto);

% bg_arr only has the first quandrant of information in it
szbg = size(bg_arr);
bg_array = zeros(szbg(1), szbg(2), 5);
for i = 1:szbg(1)
    for j = 1:szbg(2)
        bg_array(i,j,:) = cell2mat(bg_arr{i,j});
    end
end

bg_top = cat(1, flip(bg_array(2:end,:,:), 1), bg_array(:,:,:));
bg_sides = cat(2, flip(bg_top(:,2:end,:), 2), bg_top);
bg_all = cat(3, bg_sides, flip(bg_sides,3));

%imagesc(bg_all(:,:,9));

szbg = size(bg_all);
for i = 1:szbg(1)
    for j = 1:szbg(2)
        bgtemp = bg_all(i,j,:);
        
        lensmodel.TL_array(i,j) = TL(bgtemp, lensmodel.layer_thicknesses, lensmodel.TL_prototype);
        
    end
end

sgrid = lensmodel.CalcSParam(f0);
        

%pad = zeros(1000,1000);

%padded = lensmodel.PadSParams(pad, f0);

phase = angle(sgrid);
mag = abs(sgrid);

%imagesc(mag);

%plot(phase(84,:));

%TestLensHornReal(lensmodel, f0);




