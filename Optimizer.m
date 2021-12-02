classdef Optimizer
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lensmodel
        tracewidth
        optimfreqs
        element_layers
        half_element_layers
        lb
        ub
    end
    
    methods
        function obj = Optimizer(lensmodel, tracewidth, freqs)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.lensmodel = lensmodel;
            obj.tracewidth = tracewidth;
            obj.optimfreqs = freqs;
            obj.element_layers = length(lensmodel.layer_thicknesses) - 1;
            obj.half_element_layers = obj.element_layers/2;
            
            if (abs(round(obj.half_element_layers) - obj.half_element_layers) > 0.01)
                warning('Warning: optimizer cannot currently handle an odd number of metallized layers');
            end
            
            % determine upper and lower bound for element sizes. This
            % the more restrictive of the trace/space, or the limit to
            % which the elements are simulated.
            
            lb = max(min(proto.Element.sizes), obj.tracewidth/obj.lensmodel.grid_dimension);
            ub = min(max(proto.Element.sizes), 1 - obj.tracewidth/obj.lensmodel.grid_dimension);
            obj.lb = lb * ones([obj.half_element_layers, 1]);
            obj.ub = ub * ones([obj.half_element_layers, 1]);
        end
        
        function TLout = OptimTL(obj, SparamGoal, guesslist)
            sz = size(guesslist);
            guessnum = sz(2);
            TLbest = TL(guesslist(:,1), obj.lensmodel.layer_thicknesses, obj.lensmodel.TL_prototype);
            penbest = +inf;
            for i = 1:guessnum
                
                guess = guesslist(:,i);
                
                TLtemp, pentemp = obj.OptimTLsingleguess(SparamGoal, guess);
                
                if (pentemp < penbest)
                    penbest = pentemp;
                    TLbest = TLtemp;
                end
                
            end
            
            TLout = TLbest;
        end
        
        function [TLout, pen] = OptimTLsingleguess(obj, SparamGoal, initial_guess)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here      
            
                       
            layers = obj.lensmodel.layer_thicknesses;
            proto = obj.lensmodel.TL_prototype;
            optfunc = @(x) PenaltyFunc(proto.Sparam(obj.optimfreqs, [x flip(x)], layers), SparamGoal);
            options = optimset('display', 'off');

            x = fmincon(optfunc,initial_guess,[],[],[],[],obj.lb,obj.ub, [], options);
            
            TLout = TL([x flip(x)], layers, proto);
            pen = optfunc(x);
        end

       
    end
    
    methods(Static)
        function pen = PenaltyFunc(SparamIn, SparamDesired)
            S21 = SparamIn.Parameters(2,1,:);
            S21desired = SparamDesired.Parameters(2,1,:);
            
            phasepenalty = cos(angle(S21) - angle(S21desired));
            
            totalpenalty = -abs(S21) .* phasepenalty;
            
            pen = sum(totalpenalty);
        end
    end
end

