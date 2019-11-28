classdef conf
    % matlab class to handle the config for WAVE-Simulation and WAVE-Inversion
    %   members getvalue and getString read variables from config
    
    properties (Access = private)
        Variables % private cell with config variables
        Entries % private cell with config variables
    end
    
    methods
        function obj=conf(configFileName)
            % constructor which opens the config file
            if nargin>0
                 buffer = fileread( configFileName ) ; 
                 buffer = regexprep( buffer, '#.*?(\r)?\n', '$1\n' ) ;
                 Config=textscan(buffer,'%s %s','Delimiter','=');
                 obj.Variables=Config{1,1};
                 obj.Entries=Config{1,2};
            end
        end      
        function value = getValue(obj,variable)
            % get Value from config file (converts strings to num)
            index=find(strcmp(obj.Variables, variable));
            value=str2num(obj.Entries{index,1});
        end
        function string = getString(obj,variable)
            % get string from config file 
            index=find(strcmp(obj.Variables, variable));
            string=strrep(obj.Entries{index,1}, ' ', '');
        end  
    end
    
end

