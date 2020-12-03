function [FittedData, BaselineRsp, PrefRsp, PrefDir, Sigma, OppResp, Error, R2 ] = FIT_Carandini( TuningResponses, varargin )

    % set the data variables 
    if nargin > 1
        Angles = varargin{1};
        UniqueAngles = unique(Angles);
        for a = 1:length(UniqueAngles)
            TuningCurve(a) = mean( TuningResponses( Angles == UniqueAngles(a) ) );
        end
    else
        TuningCurve = TuningResponses;
        numAngles = length( TuningCurve );
        Angles = linspace( 360/numAngles, 360, numAngles );
        UniqueAngles = Angles;
    end
    
    maxBW = 90;
    minBW = 10;
    
    % guess initial settings
    [PD,~,PDx,OPPx] = TT_PreferredDirection( TuningCurve );
    BW = TT_BandWidth( TuningCurve );
    if BW > maxBW; BW = maxBW; end
    if BW < minBW; BW = minBW; end
    SG = BW/1.14;
    PR = TuningCurve(PDx) - min(TuningCurve);
    OR = TuningCurve(OPPx) - min(TuningCurve);
    BS = min(TuningCurve);
%     BW = TT_BandWidth( FittedData );
%     OSI = TT_OrientationSelectivityIndex( FittedData );
%     DI = TT_DirectionIndex( FittedData );
    % fit tuning curve
    [BaselineRsp, PrefRsp, PrefDir, Sigma, OppResp, FittedData, Error, R2] = ... 
            otfit_carandini( Angles, BS, PR, PD, 30, ...
                'widthint', [minBW maxBW], ...
                'Rpint', [mean(TuningCurve) 3*PR], ...
                'Rnint', [BS 3*PR], ...
                'spontint', [BS BS], ...
                'data', TuningResponses );
    
%     % display result
%     figure;
%     plot( Angles, TuningResponses, 'xr' );
%     hold on;
%     plot( UniqueAngles, TuningCurve, 'd:g' );
%     plot( 1:360, FittedData, 'color', [0 0 1], 'linewidth', 2 );
%     BW = TT_BandWidth( FittedData );
%     OSI = TT_OrientationSelectivityIndex( FittedData );
%     DI = TT_DirectionIndex( FittedData );
%     disp([ 'PD=' num2str(PrefDir) ', BWest=' num2str(Sigma*1.14) ...
%         ', BW=' num2str(BW) ', OSI=' num2str(OSI) ', DI=' num2str(DI) ]); 
%     disp([ 'Error = ' num2str(Error)]);
        
end