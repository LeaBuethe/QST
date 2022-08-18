

chAssign = [2,1,3];

if exist('X1','var')
        % make sure all have the same number of pulses. 
         if ~isequal(size(X1,1),size(X2,1),size(X3,1))
             X1 = X1(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
             X2 = X2(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
             X3 = X3(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
         end
 
        % set which channel ist the target channel etc
        quadratures = zeros([size(X1) 3]);
        quadratures(:,:,:,1) = X1;
        quadratures(:,:,:,2) = X2;
        quadratures(:,:,:,3) = X3;
        Xtg = quadratures(:,:,:,chAssign(1));  
        XpsFast = quadratures(:,:,:,chAssign(2));
        XpsSlow = quadratures(:,:,:,chAssign(3));
        clear('quadratures'); 
end
%%
 phaseselection = -pi:0.5:pi ;
%phaseselection = 0;


% radius and thickness for amplitude selection
r = 10;
d = 2;
w = 0.5;
%normalisation of the Quadratures from postselection channels
[nPsFast,nPsSlow,ntg] = nPhotons(O1,O2,O3);
At=sqrt(2)*10*sqrt((nPsFast+nPsSlow)*ntg)/(1+nPsFast+nPsSlow);
% nPsMean=(nPsFast+nPsSlow)/2;
% O1norm = O1*sqrt(nPsMean/nPsFast);
% O2norm = O2*sqrt(nPsMean/nPsSlow);


%[n1,n2,n3] = nPhotons(O1,O2,O3);


%nTg = nPhotons(X2,X2,X2);
[nDs,nDsStd,lengthSelect] = deal(length(phaseselection));

for i  = 1:length(phaseselection)  

       selParams = struct('Type','phase','Position',[phaseselection(i),w,r,d]);
       %selParams = struct('Type','phaseAndAmplitude','Position',[0.5,0.01,r,d]);
 [~,~,iSelect] = selectRegion(O1,O2,O3,oTheta,selParams);%,'Plot','show','Filename',['testTimes2' num2str(phaseselection(i))]);

% 
      
        selO3=O3(iSelect);
        selO2=O2(iSelect);
        selO1=O1(iSelect);
        [n1sel,n2sel,n3sel] = nPhotons(selO1,selO2,selO3);
        nDs1(i) = n1sel;
        nDs2(i) = n2sel;
        nDs3(i) = n3sel;
        nOtg = mean(nDs3);
        Aps(i)= At*(1+(n1sel+n2sel))/sqrt(2*((n1sel+n2sel)*ntg));

      %  At(i)=sqrt(2)*10*sqrt((n1sel+n2sel)*ntg)/(1+n1sel+n2sel);
         %nDsStd(i) = std(nValues);
end
%%

for i  = 1:length(phaseselection)  
     r=Aps(i);
     % r=2*Aps(i)-10;
       selParams = struct('Type','phaseAndAmplitude','Position',[phaseselection(i),w,r,d]);
       %selParams = struct('Type','phaseAndAmplitude','Position',[0.5,0.01,r,d]);
 [selX,selTheta,iSelect] = selectRegion(O1,O2,O3,oTheta,selParams);%,'Plot','show','Filename',['testTimes2' num2str(phaseselection(i))]);
 
 
    [nValues] = deal(zeros(10,1));
    for iN=1:10
        try
            uniformX = seriesUniformSampling(selX,thetaMiraSel,'NBins',100); %thetaMira oder selTheta?
        catch
            warning(['Problem with uniformSampling.', ...
                'Use X without uniform sampling.']);
            uniformX = selX;
        end  
        [nValues(iN),~,~] = nPhotons(uniformX,uniformX,uniformX);

    end 
         %[nValues,~,~] = nPhotons(selX,selX,selX);
   %lengthSelect(i)=length(iSelect);
    nDs(i) =  mean(nValues);
    nDsStd(i) = std(nValues);
end
%%
%Selection of 'fullcircle' in selectRegion for getting a mean Photonnumber
%nDs_amplitudeSelection 

selParams2 = struct('Type','fullcircle','Position',[r,d]);

      % selParams = struct('Type','phase','Position',[phaseselection(i),0.5]);
         [selX,selTheta,iSelect] = selectRegion(O1,O2,O3,oTheta,selParams);%,'Plot','show','Filename',['testTimes10' num2str(phaseselection(i)) '-assessTheta']);
%       
  %   [selX,selTheta,iSelect] = selectRegionAroundZero(O1,O2,O3,oTheta,oThetaMira,selParams2);
    %  thetaMiraSel = oThetaMira(iSelect);
    thetaMiraSel = mod(selTheta, 2*pi);
     % compute photon number of postselected quadratures in the doubleslit
    
    [nValues] = deal(zeros(10,1));
    for iN=1:10
        try
            uniformX = seriesUniformSampling(selX,thetaMiraSel,'NBins',100); %thetaMira oder selTheta?
        catch
            warning(['Problem with uniformSampling.', ...
                'Use X without uniform sampling.']);
            uniformX = selX;
        end  
        [nValues(iN),~,~] = nPhotons(uniformX,uniformX,uniformX);

    end 
    
    nDs_amplSel = mean(nValues);
    nDsStd_amplSel = std(nValues);    
    
    %%


figure(13);
% %hold on;

errorbar(phaseselection,nDs./nDs_amplSel,nDsStd,'o-','DisplayName',['-r-' num2str(r) '-d-' num2str(d) '-w-' num2str(w)]);
%errorbar(phaseselection,nDs,nDsStd,'o-','DisplayName',['-r-' num2str(r) '-d-' num2str(d) '-w-' num2str(w)]);
%errorbar(phaseselection,nDs./nDs_amplSel,nDsStd,'o-');
 %ylim([1.1 1.8]);
 xlabel('postselected phase');
 ylabel('postselected photon number');
 legend();
 grid on;
% title(['Position ',num2str(positions(j)),'mm']);
% savefig([num2str(positions(j)) 'mm-photonNumberVsPhase.fig']);
% print([num2str(positions(j)),'mm-photonNumberVsPhase.png'],'-dpng','-r300');
% clf();
        
% figure(2);
% plot(phaseselection,lengthSelect);


%fit the phaseCheckPlotting
% figure(16);
% hold on;
% %fit with a sine
% sinFkt = 'a*sin(x-b)-c';
% startPoints=[0.1 0 1];
% nDsNorm=nDs./nDs_amplSel;
% f=fit(phaseselection', nDsNorm',sinFkt,'Lower',[0,0,-3],'Upper',[10,2*pi,3], 'Start', startPoints);
% plot(f,phaseselection,nDs./nDs_amplSel,'o-');


%Plots for variating the postselection range of phase and amplitude

%title('phasewidth 0.5  amplitude r=14, d=1');
savefig('test15,500mm_withAps.fig');