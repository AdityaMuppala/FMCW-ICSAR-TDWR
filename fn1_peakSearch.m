function [fhatPeak,xPeakPos,yPeakPos] = fn1_peakSearch(fhat,x2,y2,downPad,UpPad)
%-- This code basically finds a rough estimate of a peak and then zooms in
%to find the exact peak location using fft interpolation around the peak's
%estimate

[~,pos] = max(abs(fhat(:)));
[indX,indY] = ind2sub(size(fhat),pos);

fhat2 = (fhat(indX-downPad:indX+downPad,indY-downPad:indY+downPad));
x2cut = x2(indX-downPad:indX+downPad);
y2cut = y2(indY-downPad:indY+downPad);

fhat3 = ftx(fty(fhat2));
fhat3pad = padarray(fhat3,[UpPad UpPad],0,'both');
fhat4 = iftx(ifty(fhat3pad));

x2Up = linspace(min(x2cut),max(x2cut),(2*downPad+1+2*UpPad));
y2Up = linspace(min(y2cut),max(y2cut),(2*downPad+1+2*UpPad));

[~,pos] = max(abs(fhat4(:)));
[indX2,indY2] = ind2sub(size(fhat4),pos);
fhatPeak = fhat4(indX2,indY2);
xPeakPos = x2Up(indX2);
yPeakPos = y2Up(indY2);

end

