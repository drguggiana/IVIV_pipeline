function [TuningCurve] = TT_FoldedTuningCurve( TuningCurve )
% function [OSI] = TT_FoldedTuningCurve( TuningCurve )
% 

    
    % get number of datapoints on full circle
    [a, nDataPoints] = size(TuningCurve);

    % fold (sum 180 halves) tuning curve to orientation space instead of direction space
    TuningCurve = TuningCurve( 1:round(nDataPoints/2) ) + ...
                  TuningCurve( (round(nDataPoints/2)+1):end );


end

