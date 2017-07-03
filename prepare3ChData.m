function [ X1, X2, X3 ] = prepare3ChData( filenameLO, filenameSIG )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

CALIBRATION_CH1 = 4.596047840078126e-05; % Ampere per Volt

%Norm ist the factor in the relation between the quadratures and the ladder
%operators: q = Norm*(a^{+} + a), p = Norm*i*(a^{+} - a)
%typical values are 1/sqrt(2) or 1/2.
Norm = 1/sqrt(2);

% Compute number of LO photons
dispstat('Computing number of LO photons ...','timestamp','keepthis',0);
[data8bitLO,configLO,~]=load8BitBinary(filenameLO,'dontsave');

% Compute quadratures for target quantum state
dispstat('Computing quadratures for target quantum state ...',...
    'timestamp','keepthis',0);
[data8bitSIG,configSIG,~]=load8BitBinary(filenameSIG,'dontsave');

for iCh = 1:3
    XLO = computeQuadratures(data8bitLO(:,:,iCh), configLO, CALIBRATION_CH1);

    % Calculate the variance piece-wise to compensate slow drifts (e.g. piezos)
    NLO = mean(var(XLO));

    X = computeQuadratures(data8bitSIG(:,:,iCh), configSIG, CALIBRATION_CH1);

    % Calibration of quadratures to vacuum state
    X = Norm * X / sqrt(NLO);
    
    switch iCh
        case 1
            X1 = X;
        case 2
            X2 = X;
        case 3
            X3 = X;
    end
end

end

