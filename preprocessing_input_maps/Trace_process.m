function map_value = Trace_process(in_trace,resp_type,polarity)

% %assume a sampling frequency of 10 kHz
% sf = 10000;
switch polarity
    
    case 1 %excitatory
        switch resp_type
            case 0
                %here just blank
                map_value = NaN;
            case 1
                %rectify negative
                in_trace(in_trace>0) = 0;
                %calculate the integrated current
                map_value = trapz(in_trace)./(length(in_trace));
            case 2
                %rectify negative
                in_trace(in_trace>0) = 0;
                %calculate the integrated current
                map_value = trapz(in_trace)./(length(in_trace));
            case 3
                %here just make them 0
                map_value = 0;
            case 4
                %rectify negative
                in_trace(in_trace>0) = 0;
                %calculate the integrated current
                map_value = trapz(in_trace)./(length(in_trace));
        end
    case 2 %inhibitory
        switch resp_type
            case 0
                %here just blank
                map_value = NaN;
%                 %calculate the integrated current, rectified
%                 in_trace(in_trace<0) = 0;
%                 map_value = trapz(in_trace)./(length(in_trace)/sf);
            case 1
%                 %calculate the integrated current
%                 map_value = trapz(in_trace)./(length(in_trace)/sf);
                %calculate the integrated current, rectified
                in_trace(in_trace<0) = 0;
                map_value = trapz(in_trace)./(length(in_trace));
            case 2
%                 %calculate the integrated current, rectified
%                 in_trace(in_trace<0) = 0;
%                 map_value = trapz(in_trace)./(length(in_trace)/sf);
                %calculate the integrated current, rectified
                in_trace(in_trace<0) = 0;
                map_value = trapz(in_trace)./(length(in_trace));
            case 3
                %calculate the integrated current, rectified
                in_trace(in_trace<0) = 0;
                map_value = trapz(in_trace)./(length(in_trace));
            case 4
                %calculate the integrated current, rectified
                in_trace(in_trace<0) = 0;
                map_value = trapz(in_trace)./(length(in_trace));
        end
end