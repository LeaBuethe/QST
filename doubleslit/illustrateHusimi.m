function [ ] = illustrateHusimi(varargin)

%% Validate and parse input arguments
parser = inputParser;
defaultState = 'coherent'; 
addParameter(parser,'State',defaultState);
defaultProjection = 'off';
addParameter(parser,'Projection',defaultProjection);
defaultDim = 2;
addParameter(parser,'Dim',defaultDim,@isnumeric);
defaultTheta = 45;  % Rotation angle in degrees
addParameter(parser,'Theta',defaultTheta,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[dim,projection,state,theta] = c{:};


fontsize = 24;
fontName = 'Arial';

q = -10:0.01:10;
p = q;
n = 15;
norm = 1;
%norm = 1/sqrt(2);
nTh = 1;
nCoh = 10;
%for theta = 0:10:90

if strcmp(state,'coherent')
    %45 angle --> q0 = p0
    q0 = sqrt(2*n)*norm;
    [ H ] = cohHusimi( q, p, n,'P0', q0, 'Q0', q0, 'Norm',norm);%, 'Theta',theta
end

if strcmp(state,'thermal')
    [ H ] = thermHusimi( q, p, n, 'Norm',norm );
end

if strcmp(state,'displacedthermal')
    %45 angle --> q0 = p0
    q0 = sqrt(2*nCoh)*norm;
    [ H ] = dtsHusimi( q, p, nTh, nCoh, 'P0', q0, 'Q0', q0, 'Norm',norm);%, 'Theta',theta
end

if strcmp(state,'displacedthermalPA')
    %45 angle --> q0 = p0
    q0 = sqrt(2*nCoh)*norm;
    [ H ] = dtsHusimi( q, p, nTh, nCoh, 'P0', q0, 'Q0', q0, 'Norm',norm, 'PhaseAveraged',true);%, 'Theta',theta
end

s =surf(q,p,H);
hold on;
shading flat;

g = 1:-0.01:0;
b = g;
r = ones(1,length(g));
map = [r' g' b'];
colormap(map);

colormap jet;

%xlabel('q/A','FontSize',fontsize,'FontName',fontName);
%ylabel('p/A','FontSize',fontsize,'FontName',fontName);
xlabel('q','FontSize',fontsize,'FontName',fontName);
ylabel('p','FontSize',fontsize,'FontName',fontName);
zlabel('Q(q,p)','FontSize',fontsize,'FontName',fontName);
title('Q(q,p)','FontSize',fontsize,'FontName',fontName);

if dim == 3
    set(gca,'LineWidth',2,'XColor',[0 0 0], 'YColor', [0 0 0],...
        'FontSize',22,'FontName',fontName);
    
    print(['3D-Husimi-' state],'-dpng','-r300');
    
end


if dim == 2
    view(2);
    axis image;
    set(gca,'LineWidth',3,'XColor',[0 0 0], 'YColor', [0 0 0],'Box','on',...
        'FontSize',18,'FontName',fontName,'XTick', [-10 -5 0 5 10], 'YTick', [-10 -5 0 5 10],...
        'TickLength',[0.1 0.025],'TickDir','In');
    set(gca,'XGrid','on','YGrid','on');
    
    % Plot axis 
             
    hx = annotation('arrow');
    set(hx,'parent',gca, ...
        'Position',[min(q) 0 max(q)+abs(min(q)) 0],'lineWidth',1.5);
    hy = annotation('arrow');
    set(hy,'parent',gca, ...
        'Position',[0 min(q) 0 max(q)+abs(min(q))],'lineWidth',1.5);  
    
    %projection
    if strcmp(projection,'on')
        cl = 'k';
    else
        cl = 'none'; % For some reason you must plot this anyway to see the axis arrows
    end
    z = get(s,'ZData');
    set(s,'ZData',z-10);
    projections = real(sum(H));
    plot(q,projections*max(q)/2/max(projections),'lineWidth',2,'Color',cl); %[1 0.8 0]
    hold on;
    
    print(['2D-Husimi-' state '-Proj-' projection '-Theta-' num2str(theta) '.png'],'-dpng','-r300');
      
end
%end
%print('2D-Wigner-Koh-manyThetas.png','-dpng');

end