% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


el = Element('C:\\Users\\Cassie\\Google Drive\\Dropbox\\Documents\\My Work\\Research\\metamaterial_lens\\code\\HFSS Models\\RO3003 touchstone');

sp = el.Sparams;
spm = el.Sparams_matrix;

myS = Sparam_interp(el,  20e9 * linspace(0.9, 1.1, 10), linspace(0.1, 0.9, 6));
myabcd = abcd_interp(el,  20e9 * linspace(0.9, 1.1, 10), linspace(0.1, 0.9, 6));
