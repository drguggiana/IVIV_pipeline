function [Bandwidth, BandwidthRight, BandwidthLeft] = TT_BandWidth( TuningCurve )
% function [Bandwidth] = TT_BandWidth( TuningCurve )
% 
% returns the half-maximum, half-bandwidth at the preferred direction of 
% the tuning curve. Assumes a full 360 degree measured tuning curve
%
% Input:
% - TuningCurve: Array with mean response per direction stimulus
%
% Ouput:
% - Bandwidth:      Half bandwidth at half maximum response
% - BandwidthRight: Half bandwidth at the clockwise side of the tuningcurve
% - BandwidthLeft:  Half bandwidth at the counter-clockwise side of the 
%                   tuningcurve
%
% Written by Pieter Goltstein
% Version 1.0: July 22nd, 2011
%

%     HalfMaxPercentage = 1/2;
    HalfMaxPercentage = 1/sqrt(2);

    % set the minimum of the tuning curve to zero
    if max(TuningCurve) > 0
        TuningCurve( TuningCurve<0 ) = 0;
    else
        warning('TwoPhotonToolbox:TuningCurveBelowZero', ...
            ['Tuning curve does not contain any positive values, ' ...
            'shifted values up by substracting the minimum...']);
        TuningCurve = TuningCurve - min(TuningCurve);
    end

    % get number of datapoints on full circle
    [a, nDataPoints] = size(TuningCurve);

    % calculate width of steps in tuning curve
    angleWidth = 360/length(TuningCurve);

    % find the location and value of the peak in the  data
    [a,b,peakIndex] = TT_PreferredDirection( TuningCurve );
    peakValue = TuningCurve(peakIndex);

    % calculate the indexes for the data-points left of the peak and right 
    % of the peak
    rightSteps = peakIndex:peakIndex+(nDataPoints);
    leftSteps = peakIndex:-1:peakIndex-(nDataPoints);
    rightSteps = mod(rightSteps-1,nDataPoints)+1;
    leftSteps = mod(leftSteps-1,nDataPoints)+1;

    % step through the right half of the tuning curve to find the first 
    % datapoint that is 1/sqrt(2) (70.7%) of peak response.
    for rs = 1:length(rightSteps)
        pt = TuningCurve( rightSteps(rs) );
        if pt < (peakValue .* HalfMaxPercentage)
            break;
        end
    end
    rightAngle = (rs-2)*angleWidth;

    % interpolate for precise bandwith estimate
    intAng = angleWidth - ...
        (( ((peakValue/sqrt(2))-TuningCurve(rightSteps(rs))) / ...
        (TuningCurve(rightSteps(rs-1))-TuningCurve(rightSteps(rs))) ) ...
        * angleWidth);

    % exact estimate of right angle:
    rightAngle = rightAngle + intAng;

    % step through the left half of the tuning curve to find the first 
    % datapoint that is 1/sqrt(2) (70.7%) of peak response.
    for ls = 1:length(leftSteps)
        pt = TuningCurve( leftSteps(ls) );
        if pt < (peakValue .* HalfMaxPercentage)
            break;
        end
    end
    leftAngle = (ls-2)*angleWidth;

    % interpolate for precise bandwith estimate
    intAng = angleWidth - ...
        (( ((peakValue.*HalfMaxPercentage)-TuningCurve(leftSteps(ls))) / ...
        (TuningCurve(leftSteps(ls-1))-TuningCurve(leftSteps(ls))) ) ...
        * angleWidth);

    % exact estimate of right angle:
    leftAngle = leftAngle + intAng;
    
    % bandwidth is defined as one-half of the difference between these two angles
    Bandwidth = (rightAngle + leftAngle) / 2;
    BandwidthRight = rightAngle;
    BandwidthLeft = leftAngle;

    % if bandwidth is bigger than 180 degrees, give max as answer
    if Bandwidth > 180
        Bandwidth = 180;
    end
    if BandwidthRight > 180
        Bandwidth = 180;
    end
    if BandwidthLeft > 180
        Bandwidth = 180;
    end
    
    if Bandwidth <= 0
        Bandwidth = NaN;
        BandwidthLeft = NaN;
        BandwidthRight = NaN;
    end
    
end

