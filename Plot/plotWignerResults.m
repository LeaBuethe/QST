function plotWignerResults(listOfParams,varargin)
%PLOTSERIESPOSTSELECTIONS

%% Validate and parse input arguments
p = inputParser;
defaultXUnit = 'ps';
addParameter(p,'XUnit',defaultXUnit,@isstr);
defaultFitType = 'gauss';
addParameter(p,'FitType',defaultFitType,@isstr);
defaultVaryAPS = false;
addParameter(p,'VaryAPS',defaultVaryAPS,@islogical);
defaultRemoveModulation = false;
addParameter(p,'RemoveModulation',defaultRemoveModulation,@islogical);
defaultRange = 0.3;
addParameter(p,'Range',defaultRange,@isvector);
defaultZeroDelay = 0;
addParameter(p,'ZeroDelay',defaultZeroDelay,@isnumeric); % in mm
defaultSmooth = false;
addParameter(p,'Smooth',defaultSmooth,@islogical);
defaultPlotrelative = false; %plots the data relative to the target photon number
addParameter(p,'Plotrelative',defaultPlotrelative,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[fitType,plotrelative,range,remMod,smooth,varyAPS,xUnit,zeroDelay] = c{:};



%% Create folder 'figures-fig'
if ~exist([pwd 'Wigner-figures-fig'],'dir')
    mkdir('Wigner-figures-fig')
end

% Constants
figurepath = 'Wigner-figures-fig/';

load('Photonnumbers.mat','nPsFast','nPsSlow','nTg'); %loads the photon numbers of each channel without postselection,...
%which was created with plotSeriesPostselections
%% Gather data
[delay,Yr,Yt,Qs,Ps,varQs,varPs,discN,meanPh,meanR,varPh,varR,meanAbsPh] = deal([]);
%[delay,Yr,Yt,Qs,Ps,varQs,varPs,discN,meanPh,meanR,varPh,varR,meanAbsPh,g1WA,g1,g1Binned] = deal([]);
for iParams = 1:length(listOfParams)
    selParams = listOfParams(iParams);
    selStr = selParamsToStr(selParams);
    foldername = ['Wignerplots-',selStr,'-remMod-',...
            num2str(remMod),'-range-',num2str(range),'-varyAPS-',num2str(varyAPS)];
    load([foldername '\Wignerresults-smooth-' num2str(smooth) '.mat'],'Delay','Q','P','varQ','varP','n',...
        'meanPhases','meanAmps','varPhases','varAmps','meanAbsPhases','g1WithAmps','g1WithoutAmps','g1WithoutAmpBinneds','g1WithAmpNorms');
    delay(iParams,:) = Delay; 
    H = length(delay);
    Radii = ones(H,1) * selParams.Position(1);
    Thicknesses = ones(H,1) * selParams.Position(2);
    Yr(iParams,:) = Radii;
    Yt(iParams,:) = Thicknesses;
    Qs(iParams,:) = Q;
    Ps(iParams,:) = P;
    varQs(iParams,:) = varQ;
    varPs(iParams,:) = varP;
    discN(iParams,:) = n;
    meanPh(iParams,:) = meanPhases;
    meanR(iParams,:) = meanAmps; 
    varPh(iParams,:) = varPhases;
    varR(iParams,:) = varAmps;
    meanAbsPh(iParams,:) = meanAbsPhases;
    g1WA(iParams,:) = g1WithAmps;
    g1(iParams,:) = g1WithoutAmps;
    g1Binned(iParams,:) = g1WithoutAmpBinneds;
    g1WAN(iParams,:) = g1WithAmpNorms;
end
[~,I] = sort(delay(:,1)); % Sort for Radii
delay = delay(I,:); %delays
Yr = Yr(I,:); %Radii
Yt = Yt(I,:); %thicknesses
Qs= Qs(I,:);
Ps= Ps(I,:);
varQs = varQs(I,:);
varPs = varPs(I,:);
discN = discN(I,:);
meanPh = meanPh(I,:);
meanR = meanR(I,:);varPh = varPh(I,:);varR = varR(I,:);meanAbsPh = meanAbsPh(I,:);
g1WA=g1WA(I,:); g1=g1(I,:); g1Binned = g1Binned(I,:); g1WAN = g1WAN(I,:);
[fitTau,fitPeak,tauError] = deal(zeros(length(I),1));
if plotrelative
    meanR = meanR./sqrt(nTg);
    Qs = Qs./sqrt(nTg);
    Ps = Ps./sqrt(nTg);
    discN = discN./nTg;
    varQs = varQs./nTg;
    varPs = varPs./nTg;
end

if varyAPS
    sel = 'k'; %the selection radius or variable
else
    sel = 'r';
end

%% Plot
%typestrVector = {'R','Q','P','RLog','meanPh','meanAbsPh','discN','varR','varPh','varQ','varP','g1','g1WA','g1Binned','g1WAN','g1toOne'};
%typestrVector = {'R','discN','varR'};
%typestrVector = {'R'};
typestrVector = {'g1','g1WA','g1Binned','g1WAN','g1toOne'};
for typeI = 1:length(typestrVector)
    plotStuff(cell2mat(typestrVector(typeI)))
end

    function [] = plotStuff(typestr)
%% Create figure
fig = figure;
filename = [figurepath,typestr,'-',fitType,'-remMod-',...
            num2str(remMod),'-range-',num2str(min(range)),'-',num2str(max(range)),'-varyAPS-',num2str(varyAPS),'-smooth-',num2str(smooth),'-plotrelative-' num2str(plotrelative),'.fig'];
%formatFigA5(fig);
figure(1);
ax = gca;
for i = 1:length(I)
    switch typestr
        case 'R'
            ys = meanR(i,:);
            if plotrelative
                ylabel(ax,'$mean amplitude <r>/\sqrt{n_{Tg}}$','Interpreter','latex'); 
            else            
                ylabel(ax,'mean amplitude <r>'); 
            end
        case 'RLog'
            ax.YScale = 'log';
            ys = meanR(i,:);
            if plotrelative
                ylabel(ax,'$mean amplitude <r>/\sqrt{n_{Tg}}$','Interpreter','latex'); 
            else            
                ylabel(ax,'mean amplitude <r>'); 
            end
        case 'Q'
            ys = Qs(i,:);
            if plotrelative
                ylabel(ax,'$<Q>/\sqrt{n_{Tg}}$','Interpreter','latex'); 
            else            
                ylabel(ax,'<Q>'); 
            end       
        case 'P'
            ys = Ps(i,:);
            if plotrelative
                ylabel(ax,'$<P>/\sqrt{n_{Tg}}$','Interpreter','latex'); 
            else            
                ylabel(ax,'<P>'); 
            end        
        case 'meanPh'
            ys = meanPh(i,:);
            ylabel(ax,'<\phi> (rad)'); 
        case 'meanAbsPh'
             ys = meanAbsPh(i,:);          
            ylabel(ax,'<|\phi|> (rad)'); 
        case 'varPh'
            ys = varPh(i,:); 
            ylabel(ax,'Variance(\phi) (rad)');
        case 'varR'    
            ys = varR(i,:);
            ylabel('Variance(r)');
        case 'varQ'
            ys = varQs(i,:);
             if plotrelative
                ylabel(ax,'$Var(Q)/n_{Tg}$','Interpreter','latex'); 
            else            
                ylabel(ax,'Var(Q)'); 
             end
         case 'varP'
             ys = varPs(i,:);
             if plotrelative
                ylabel(ax,'$Var(P)/n_{Tg}$','Interpreter','latex'); 
            else            
                ylabel(ax,'Var(P)'); 
             end       
        case 'discN'
            ys = discN(i,:);
            if plotrelative
                ylabel(ax,'$n/n_{Tg}$','Interpreter','latex'); 
            else            
                ylabel(ax,'n'); 
            end
        case 'g1'
            ys = g1(i,:);
            ylabel(ax,'g^{(1)} only phase'); 
        case 'g1toOne'
            ys = g1(i,:)/max(g1(i,:));
            ylabel(ax,'g^{(1)} only phase');
        case 'g1WA'
            ys = g1WA(i,:);
            ylabel(ax,'g^{(1)} with amplitude'); 
        case 'g1Binned'
            ys = g1Binned(i,:);
            ylabel(ax,'g^{(1)} only phase, binned'); 
        case 'g1WAN'
            ys = g1WAN(i,:);
            ylabel(ax,'g^{(1)} with amplitude, norm.'); 
    end %typestr     
        
    plot(ax,delay(i,:),ys,'o-','DisplayName',[sel ' = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
    hold on;
    x = delay(i,:);
    if ~any(~isnan(ys)) %if there are only nans
        ys = zeros(1,H);
    end
            ys = transpose(csaps(x,ys,0.00001,x));  
            ys = ys';
    switch fitType
        case 'gauss'
            g = fittype(@(a,b,c,d,x) a*exp(-pi/2*((x-b)/c).^2)+d);                   
            b0 = zeroDelay;
            c0 = 0.5*max(x);
            d0 = mean(ys(end-5:end)); % saturation value
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
            [res,gof,~] = fit(x',ys',g,'StartPoint',[a0 b0 c0 d0]); 
            fitTau(i) = res.c;
            level = 2*tcdf(-1,gof.dfe);
            m = confint(res,1-level); 
            tauError(i) = m(end,3) - res.c;
            fitPeak(i) = res(res.b);  
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),res(min(delay(i,:)):1:max(delay(i,:))),'r','DisplayName','');
        case 'exponential'
            ex = fittype(@(m,c,x) m - x/c);                   
            b0 = zeroDelay;
            c0 = 0.5*max(x);                
            d= mean(ys(end-5:end));
            ys = ys-d; 
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay)));
            sig = sign(a0);
            ys = abs(ys);
            a0 = abs(a0);
            m0 = log(a0) + b0/c0;
            [res,gof,~] = fit(abs(x)',log(ys)',ex,'StartPoint',[m0 c0]); 
            fitTau(i) = res.c;
            level = 2*tcdf(-1,gof.dfe);
            m = confint(res,1-level); 
            tauError(i) = m(end,2) - res.c;         
            fitPeak(i) = sig*exp(res(b0))+d;  
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),sig*exp(res(abs(min(delay(i,:)):1:max(delay(i,:)))))+d,'r','DisplayName',''); 
        case 'exp2'
            g = fittype(@(a,b,c,d,x) a*exp(-((x-b)/c)) +d);                   
            b0 = zeroDelay;
            c0 = 0.5*max(x);
            d0 = mean(ys(end-5:end)); % saturation value
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
            [res,gof,~] = fit(abs(x)',ys',g,'StartPoint',[a0 b0 c0 d0]); 
            fitTau(i) = res.c;
            level = 2*tcdf(-1,gof.dfe);
            m = confint(res,1-level); 
            tauError(i) = m(end,3) - res.c;
            fitPeak(i) = res(res.b);  
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),res(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');
        case 'sech'
            se = fittype(@(a,b,c,d,x) a./(cosh(-pi/2*(x-b)/c)).^2+d);                   
            b0 = zeroDelay;
            c0 = 0.5*max(x);
            d0 = mean(ys(end-5:end)); % saturation value
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
            [res,gof,~] = fit(x',ys',se,'StartPoint',[a0 b0 c0 d0]); 
            fitTau(i) = res.c;
            level = 2*tcdf(-1,gof.dfe);
            m = confint(res,1-level); 
            tauError(i) = m(end,3) - res.c;
            fitPeak(i) = res(res.b);  
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),res(min(delay(i,:)):1:max(delay(i,:))),'r','DisplayName','');
        case 'power-law'
            g = fittype(@(a,b,alpha,x) a*abs(2*(x-b)).^(-alpha));                   
            b0 = zeroDelay;
            alpha0 = 0.22;
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay)));
            [res,gof,~] = fit(x',ys',g,'StartPoint',[a0 b0 alpha0]); 
            fitTau(i) = res.alpha;
            [se]= getStandardErrorsFromFit(res,gof,'method1');   
            tauError(i) = se(3);
            fitPeak(i) = res(res.a);  
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),res(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');           
         case 'lorentz'
             % starting values 
            p02 = zeroDelay;
            p03 = 0.5*max(x);
            C = mean(ys(end-5:end)); % saturation value 
            peakHeight = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))); 
            p01 = (peakHeight- C)*p03;
            p0 = [p01 p02 p03 C];
             [res,params,RESNORM, RESIDUAL, JACOBIAN] = lorentzfit(x',ys',p0);
            p1 = params(1);  p2 = params(2);  p3 = params(3); C = params(4);
            fitTau(i) = 2*sqrt(p3);  %FWHM
            fitPeak(i) = p1/p3 + C;  
            plot(ax,delay(i,:),res,'r','DisplayName','');                   
        case 'voigt' 
%                     fo = fitoptions('Method','NonlinearLeastSquares', ...
%                             'StartPoint',[50, 5]); %gamma, sigma
%                     zeroDelay = 200;  %set this...  
%                     ysFit = ys/max(ys);
%                     ft = fittype( @(gam,sig,peakVal,x) myvoigt(x, peakVal, gam, sig) ,'problem','peakVal','options',fo);    
%                     [res,gof,~] = fit(x',ysFit,ft,'problem',zeroDelay);
%                     resPlot= res(min(delay(i,:)):1:max(delay(i,:)));
%                     resPlot = resPlot * max(ys);
%                     plot(ax,min(delay(i,:)):1:max(delay(i,:)),resPlot,'r','DisplayName','');              
            initGuess1 = [0,30, 1]; %peak x, gamma, sigma
            [estimates1, model1] = voigtfit(x, ys, initGuess1, [min(x), max(x)]);
            gamma = estimates1(2);
            sigma = estimates1(3);
            res = myvoigt(x, 0, gamma, sigma );%estimates1(1)
            res = res';
            c = abs(res\ys');
            res = res*c;             
            plot(ax,x,res,'r','DisplayName','');
            FWHM_Gauss = 2*sigma*sqrt(2*log(2));
            FWHM_Lorentz = 2*gamma;
            fitTau(i) = 0.5346*FWHM_Lorentz + sqrt(0.2166*FWHM_Lorentz^2 + FWHM_Gauss^2);
            fitPeak(i) = max(res);
            % https://en.wikipedia.org/wiki/Voigt_profile
            % Olivero, J. J.; R. L. Longbothum (February 1977). 
%                     "Empirical fits to the Voigt line width: A brief review". 
%                     Journal of Quantitative Spectroscopy and Radiative Transfer. 
%                     17 (2): 233???236. Bibcode:1977JQSRT..17..233O. doi:10.1016/0022-4073(77)90161-3. ISSN 0022-4073.
        case 'noFit'
            fitPeak(i) = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay)));
    end %fitType
end %iloop
hold off;
if strcmp(fitType,'noFit')
    legend('location','best');
else
    f=get(ax,'Children');
    index = length(f)-((1:length(I))-1).*2;
    legend(f(index),'location','northeast');
end
xlabel(ax,['Delay (' xUnit ')']);
fig = figure(1);              
set(fig,'Color','w','Units','centimeters','Position',[1,1,45,30],'PaperPositionMode','auto');
graphicsSettings;
set(ax,'FontSize',38);
%% Write figure to file and close it
if ~isempty(filename)
    savefig(fig,filename);
    print([filename '.png'],'-dpng');
    close all;
end

% %% plot fit results vs ring radius and thickness
if any(fitTau)
    YrPlot = Yr(:,1);
    [YrPlot,Ir]= sort(YrPlot);
    fitTau = real(fitTau(Ir));
    tauError = tauError(Ir);
    errorbar(YrPlot,fitTau,tauError,'o-','Linewidth',2);
    xlabel('A_{ps} set for postselection');
    switch fitType
        case 'gauss'
            ylabel(['\tau_c of Gaussian (' xUnit ')']);
        case 'voigt'
            ylabel(['FWHM of Voigt Profile (' xUnit ')']);
        case 'exponential'
            ylabel(['\tau_c of exp. function (' xUnit ')']);
        case 'exp2'
            ylabel(['\tau_c of exp. function (' xUnit ')']);
        case 'sech'
            ylabel(['\tau_c of sech function (' xUnit ')']);
         case 'lorentz'
            ylabel(['FWHM of Lorentz Profile (' xUnit ')']);
         case 'power-law'
            ylabel('\alpha of power-law Profile');
    end
    graphicsSettings;
    ax = gca;
    set(ax,'FontSize',30);
    title(typestr);
    savefig([filename '-fitWidthsVsRadius.fig']);
    print([filename '-fitWidthsVsRadius.png'],'-dpng');
    close all;
end

if any(fitPeak)
    YrPlot = Yr(:,1);
    [YrPlot,Ir]= sort(YrPlot);
    fitPeak = real(fitPeak(Ir));
    plot(YrPlot,fitPeak,'o-');
    xlabel('A_{ps} set for postselection');
    ylabel('Peakheight');
    graphicsSettings;
    ax = gca;
    set(ax,'FontSize',30);
    title(typestr);
    savefig([filename '-PeaksVsRadius.fig']);
    print([filename '-PeaksVsRadius.png'],'-dpng');
    close all;
end
% 
% if any(fitTau)
%     YtPlot = Yt(:,1);
%     [YtPlot,It]= sort(YrPlot);
%     fitTau = fitTau(It);
%     tauError = tauError(It);
%     errorbar(YtPlot,fitTau,tauError,'o-','Linewidth',2);
%     xlabel('Ring Thickness set for postselection');
%     switch fitType
%         case 'gauss'
%             ylabel(['\tau_c of Gaussian (' xUnit ')']);
%         case 'voigt'
%             ylabel(['FWHM of Voigt Profile (' xUnit ')']);
%          case 'exponential'
%             ylabel(['\tau of exp. function (' xUnit ')']);
%     end
%     graphicsSettings;
%     title([filename '-' typestr '-' fitType '-fitWidthsVsThickness']);
%     savefig([filename '-' typestr '-' fitType '-fitWidthsVsThickness.fig']);
%     print([filename '-' typestr '-' fitType '-fitWidthsVsThickness.png'],'-dpng');
%     close all;
% end
% 
% if any(fitPeak)
%     fitPeak = fitPeak(It);
%     plot(YtPlot,fitPeak,'o-');
%     xlabel('Ring Thickness set for postselection');
%     ylabel('Peakheight');
%     graphicsSettings;
%     title([filename '-' typestr '-' fitType '-PeaksVsThickness']);
%     savefig([filename '-' typestr '-' fitType '-PeaksVsThickness.fig']);
%     print([filename '-' typestr '-' fitType '-PeaksVsThickness.png'],'-dpng');
%     close all;
% end
% plot(Yt(:,1),fitWidth,'o-');
% plot(Yt(:,1),fitPeak,'o-');
    end %plotStuff

end

