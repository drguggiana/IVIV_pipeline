function [PreferredDirection, OppositeDirection, PreferredIndex, ...
                    OppositeIndex] = TT_PreferredDirection( TuningCurve )
% function [PreferredDirection, OppositeDirection, PreferredIndex, ...
%                     OppositeIndex] = TT_PreferredDirection( TuningCurve )
%
% Returns the direction that gives rise to the largest response and the
% opposite of this direction. Assumes a full 360 degree measured
% tuningcurve.
%
% Input:
% - TuningCurve: Array with mean response per direction stimulus
%
% Ouput:
% - PreferredDirection: Direction with largest response
% - OppositeDirection:  Opposite direction of preferred direction
% - PreferredIndex:     Index of the tuningcurve array at the preferred dir
% - OppositeIndex:      Index of the tuningcurve array at the opposite dir
%
% Written by Pieter Goltstein
% Version 1.0: July 22nd, 2011
%
    
    % get number of datapoints on full circle
    [a, nDataPoints] = size(TuningCurve);
    Directions = (1:nDataPoints) * (360/nDataPoints);

    % find array indices with the largest response
    PreferredIndex = find( TuningCurve == max(TuningCurve) );
    
    % if more than one peak has been found, take the mean of the
    % neighboring data points in account
    if length(PreferredIndex) > 1
        
        max3peak = zeros(length(PreferredIndex),1);
        for d = 1:length(PreferredIndex)
            Neighbors = mod( ((PreferredIndex(d)-1):1:(PreferredIndex(d)+1))-1, length(TuningCurve) ) + 1;
            max3peak(d) = mean( TuningCurve( Neighbors ) );
        end
        MaxMaxIndx = find( max3peak == max(max3peak) );
        
        % if this doesnt decide, do it with 2 neighbors on each side
        if length(MaxMaxIndx) > 1
            
            max5peak = zeros(length(PreferredIndex),1);
            for d = 1:length(PreferredIndex)
                Neighbors = mod( ((PreferredIndex(d)-2):1:(PreferredIndex(d)+2))-1, length(TuningCurve) ) + 1;
                max5peak(d) = mean( TuningCurve( Neighbors ) );
            end
            MaxMaxIndx = find( max5peak == max(max5peak) );

            % if this doesnt decide, spit out a warning and make a random
            % choice between the maximum responses
            if length(MaxMaxIndx) > 1
                warning('TwoPhotonToolbox:MultiplePreferredDirections', ...
                    'Found multiple identical peaks, cannot decide preferred direction, taking random pick...')
                MaxMaxIndx = ceil(rand(1)*length(PreferredIndex));
                if MaxMaxIndx == 0
                    MaxMaxIndx = 1;
                end
            end                
        end
        
        PreferredIndex = PreferredIndex( MaxMaxIndx );
    end
    
    OppositeIndex = mod( round(PreferredIndex+(length(TuningCurve)/2))-1, length(TuningCurve) )+1;
    
    PreferredDirection = Directions( PreferredIndex );
    
    OppositeDirection = Directions( OppositeIndex );
    
end