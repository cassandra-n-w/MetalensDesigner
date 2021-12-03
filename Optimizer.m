classdef Optimizer < handle
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
        TLlist
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
            
            obj.TLlist = TL.empty(0);
            
            if (abs(round(obj.half_element_layers) - obj.half_element_layers) > 0.01)
                warning('Warning: optimizer cannot currently handle an odd number of metallized layers');
            end
            
            % determine upper and lower bound for element sizes. This
            % the more restrictive of the trace/space, or the limit to
            % which the elements are simulated.
            proto = lensmodel.TL_prototype;
            
            lb = max(min(proto.Element.sizes), obj.tracewidth/obj.lensmodel.grid_dimension);
            ub = min(max(proto.Element.sizes), 1 - obj.tracewidth/obj.lensmodel.grid_dimension);
            obj.lb = lb * ones([obj.half_element_layers, 1]);
            obj.ub = ub * ones([obj.half_element_layers, 1]);
        end
        
        function Optim360(obj)
            % ONLY FOR SINGLE FREQUENCY OPTIMIZATION
            ang = linspace(0,360,31);
            
            for i = 1:length(ang)
                
                
                goal = exp(1.0i * ang(i)*pi/180);
                
                TLtemp = obj.OptimTL(goal, obj.GenGuesses());
                obj.TLlist(i) = TLtemp;
            end
        end
        
        function TLout = OptimTL(obj, S21Goal, guesslist)
            sz = size(guesslist);
            guessnum = sz(1);
            TLbest = TL(guesslist(:,1), obj.lensmodel.layer_thicknesses, obj.lensmodel.TL_prototype);
            penbest = +inf;
            for i = 1:guessnum
                
                guess = guesslist(i,:);
                
                [TLtemp, pentemp] = obj.OptimTLsingleguess(S21Goal, guess);
                
                if (pentemp < penbest)
                    penbest = pentemp;
                    TLbest = TLtemp;
                end
                
            end
            
            penbest
            
            TLout = TLbest;
        end
        
        function [TLout, pen] = OptimTLsingleguess(obj, S21Goal, initial_guess)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here      
            
                       
            layers = obj.lensmodel.layer_thicknesses;
            proto = obj.lensmodel.TL_prototype;
            optfunc = @(x) Optimizer.PenaltyFunc(proto.SParam(obj.optimfreqs, [x flip(x)], layers), S21Goal);
            options = optimset('display', 'off');

            x = fmincon(optfunc,initial_guess,[],[],[],[],obj.lb,obj.ub, [], options);
            
            TLout = TL([x flip(x)], layers, proto);
            pen = optfunc(x);
        end

       
        function guesses = GenGuesses(obj)
            guessarr = [.15, .4, .85];
            
            
            guesses = [[.1 .1 .1 .1 .1] ; [.9 .9 .9 .9 .9] ; [.1 .9 .1 .9 .1]; [.9 .1 .9 .1 .9]];
        end
    end
    
    methods(Static)
        function pen = PenaltyFunc(SparamIn, S21Desired)
            S21 = SparamIn(2,1,:);
            
            phasepenalty = cos(angle(S21) - angle(S21Desired));
            
            totalpenalty = -abs(S21) .* phasepenalty;
            
            pen = mean(totalpenalty);
        end
        
    end
end

