classdef SemiautoMeas
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Z
        X
        Y
        Sparam
        Frequencies
    end

    methods
        function obj = SemiautoMeas(path, file, Xspacing)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            filedata = matfile(fullfile(path, file));

            newstr = split(file, "_");
    
            newstr = regexprep(newstr,'[^0-9]','');

            %obj.Z = str2num(newstr{1}) + str2num(newstr{2})/100;
            %obj.X = str2num(newstr{3}) + str2num(newstr{4})/100;
            obj.Y = str2num(newstr{1}) + str2num(newstr{2})/100;

            obj.Sparam = filedata.data;
            obj.Frequencies = filedata.freq * 1e9;
            x_count = size(obj.Sparam);
            x_count = x_count(1);

            obj.X = 0:Xspacing:((x_count-1)*Xspacing);
            obj.Z = 0;

            obj.Sparam = reshape(obj.Sparam, x_count, 1, []);
        end

        function outputArg = select_data(obj, freqs)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            idx = interp1(obj.Frequencies, 1:length(obj.Frequencies),freqs,'nearest'); 
            outputArg = obj.Sparam(:,:,round(idx));
        end
    end
end