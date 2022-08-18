function[thetaCorr,thetaMiraCorr,XtgCorr,XpsSlowCorr,XpsFastCorr]=corrBoundariesOfPiezosegments(theta,thetaMira,Xtg,XpsSlow,XpsFast);
[seg,piezoSeg]=size(theta);
thetaCorr=theta;
thetaCorr(102952:seg,:)=[];
thetaCorr(1:20800,:)=[];

thetaMiraCorr=thetaMira;
thetaMiraCorr(102952:seg,:)=[];
thetaMiraCorr(1:20800,:)=[];
%%
Xtg = reshape(Xtg,size(theta));
XtgCorr=Xtg;
XtgCorr(102952:seg,:)=[];
XtgCorr(1:20800,:)=[];

XpsSlow = reshape(XpsSlow,size(theta));
XpsSlowCorr=XpsSlow;
XpsSlowCorr(102952:seg,:)=[];
XpsSlowCorr(1:20800,:)=[];

XpsFast = reshape(XpsFast,size(theta));
XpsFastCorr=XpsFast;
XpsFastCorr(102952:seg,:)=[];
XpsFastCorr(1:20800,:)=[];
end