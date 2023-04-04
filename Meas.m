classdef Meas
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
        function obj = Meas(path, file)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            filedata = load(fullfile(path, file));

            newstr = split(file, "_");
    
            newstr = regexprep(newstr,'[^0-9]','');

            obj.Z = str2num(newstr{1}) + str2num(newstr{2})/100;
            obj.X = str2num(newstr{3}) + str2num(newstr{4})/100;
            obj.Y = str2num(newstr{5}) + str2num(newstr{6})/100;
            obj.Sparam = filedata(:,2);
            obj.Frequencies = filedata(:,1);
        end

        function outputArg = select_data(obj, freqs)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            idx = interp1(obj.Frequencies, 1:length(obj.Frequencies),freqs,'nearest'); 
            outputArg = obj.Sparam(idx);
        end
    end
end