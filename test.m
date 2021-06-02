el = Element('C:\\Users\\Cassie\\Google Drive\\Dropbox\\Documents\\My Work\\Research\\metamaterial_lens\\code\\HFSS Models\\RO3003 touchstone');

sp = el.Sparams;
spm = el.Sparams_matrix;

myS = Sparam_interp(el,  20e9 * linspace(0.9, 1.1, 10), linspace(0.1, 0.9, 6));
myabcd = abcd_interp(el,  20e9 * linspace(0.9, 1.1, 10), linspace(0.1, 0.9, 6));
