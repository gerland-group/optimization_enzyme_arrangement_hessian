load '../E_prof.out';
load '../returns.out';
load '../flux_comp.out';

fname = '../par.out';
fid = fopen(fname);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);

pars = jsondecode(str);

Nx=str2double(pars.Nx);
Ny=str2double(pars.Ny);

sources = fieldnames(pars.sources);
if(isfield(pars,'absorbers'))
    absorbers = fieldnames(pars.absorbers);
end
%definig the frame
m=size(E_prof(:,1),1);
n=size(E_prof(1,:),2);

Elast=E_prof((m-Ny+1):m,1:n);

%axis tight;
set(0,'DefaultAxesFontSize',6)
set(0,'DefaultAxesFontName', 'Helvetica')
set(0,'DefaultTextFontName', 'Helvetica')

%image
%figure('Units','inches','Position',[0,20,4,4],'PaperUnits', 'Inches', 'PaperSize', [5, 5],'PaperPosition',[0 0 5 5]);
figure('Units','inches','Position',[20,10,6.5,6.5],'PaperUnits', 'Inches', 'PaperSize', [6.5 , 4.5],'PaperPosition',[0 0 6.5 6.5]);
colormap('bone');
colormap(flipud(colormap));

%MACROS FOR STYLE
%line_width=0.7;
line_width=0.9;
x_com_size=0.4;
y_com_size=0.4;
line_style_sources='--';
%positions to compare [row column]
pos_1=[10 9];
color_1 = 'b';
pos_2=[15 24];
color_2 = 'r';
pos_3=[30 14];
color_3 = [0 0.5 0];

e_pos_1_last = Elast(pos_1(1),pos_1(2));
e_pos_2_last = Elast(pos_2(1),pos_2(2));
e_pos_3_last = Elast(pos_3(1),pos_3(2));
e_pos_max = max([e_pos_1_last,  e_pos_2_last, e_pos_3_last]);
E_T_last = sum(Elast(:))/(Nx*Ny);
J_comp_uniform_max = max((flux_comp(:,2)./flux_comp(:,3)));
J_comp_uniform_min = min((flux_comp(:,2)./flux_comp(:,3)));
J_comp_cluster_max = max((flux_comp(:,2)./flux_comp(:,4)));
J_comp_cluster_min = min((flux_comp(:,2)./flux_comp(:,4)));

%FRAMES
%FRAME WITH FEW Enzymes
ax1 = subplot(3,3,1);
i=floor(m/(Ny*12));
E_frame=E_prof((i-1)*Ny+1:i*Ny,1:n);
imagesc(E_frame);
xlabel('x');
xticks([0.5 Nx+0.5]);
xticklabels({0,Nx});
set (gca,'Ydir','normal');
ylabel('y');
yticks([0.5 Ny+0.5]);
yticklabels({0,Ny});
set(ax1, 'YLim', [0.5 Ny+0.5],'XLim', [0.5 Nx+0.5]);

E_T = sum(E_frame(:));

for j = 1:numel(sources)
    sourcej=pars.sources.(sources{j});
    A_0=str2double(sourcej.A_0);
    SPx=str2double(sourcej.x_pos);
    SPy=str2double(sourcej.y_pos);
    SSx=str2double(sourcej.x_size);
    SSy=str2double(sourcej.y_size);
    rectangle('Position',[SPx+0.5 SPy+0.5 SSx SSy],'Linewidth',line_width,'LineStyle',line_style_sources);
end

if(isfield(pars,'absorbers'))
        for j = 1:numel(absorbers)
            absorberj=pars.absorbers.(absorbers{j});
            APx=str2double(absorberj.x_pos);
            APy=str2double(absorberj.y_pos);
            ASx=str2double(absorberj.x_size);
            ASy=str2double(absorberj.y_size);
            rectangle('Position',[APx+0.5 APy+0.5 ASx ASy],'Linewidth',line_width);
        end
end

if E_T>0
    [y_grid,x_grid] = ndgrid(1:size(E_frame,1),1:size(E_frame,2));
    y_com = sum(y_grid(:).*E_frame(:))/E_T;
    x_com = sum(x_grid(:).*E_frame(:))/E_T;
    
    %rectangle('Position',[x_com-x_com_size/2 y_com-y_com_size/2 x_com_size y_com_size],'Curvature',[1 1],'Linewidth',line_width, 'FaceColor','w');
end

rectangle('Position',[0.5 0.5 Nx Ny],'Linewidth',line_width);
rectangle('Position',[pos_1(2)-0.5 pos_1(1)-0.5 1 1],'LineWidth',line_width,'EdgeColor',color_1);
rectangle('Position',[pos_2(2)-0.5 pos_2(1)-0.5 1 1],'LineWidth',line_width,'EdgeColor',color_2);
rectangle('Position',[pos_3(2)-0.5 pos_3(1)-0.5 1 1],'LineWidth',line_width,'EdgeColor',color_3);

caxis([0 max(max(Elast))]);
t1=title(['E_T=',num2str(floor(E_T/(Nx*Ny)))],'FontWeight','normal');
%SECOND FRAME
%ax2 = subplot(2,3,2);
ax2 = subplot(3,3,2);
E_frame=Elast;
imagesc(E_frame);
xlabel('x');
xticks([0.5 Nx+0.5]);
xticklabels({0,Nx});
set (gca,'Ydir','normal');
ylabel('y');
yticks([0.5 Ny+0.5]);
yticklabels({0,Ny});
set(ax2, 'YLim', [0.5 Ny+0.5],'XLim', [0.5 Nx+0.5]);

E_T = sum(E_frame(:));

for j = 1:numel(sources)
    sourcej=pars.sources.(sources{j});
    A_0=str2double(sourcej.A_0);
    SPx=str2double(sourcej.x_pos);
    SPy=str2double(sourcej.y_pos);
    SSx=str2double(sourcej.x_size);
    SSy=str2double(sourcej.y_size);
    rectangle('Position',[SPx+0.5 SPy+0.5 SSx SSy],'Linewidth',line_width,'LineStyle',line_style_sources);
end

% Enzyme density center of mass
% if E_T>0
%     [y_grid,x_grid] = ndgrid(1:size(E_frame,1),1:size(E_frame,2));
%     y_com = sum(y_grid(:).*E_frame(:))/E_T;
%     x_com = sum(x_grid(:).*E_frame(:))/E_T;
%         
%     rectangle('Position',[x_com-x_com_size/2 y_com-y_com_size/2 x_com_size y_com_size],'Curvature',[1 1],'Linewidth',line_width, 'FaceColor','w');
% end

rectangle('Position',[0.5 0.5 Nx Ny],'Linewidth',line_width);
rectangle('Position',[pos_1(2)-0.5 pos_1(1)-0.5 1 1],'LineWidth',line_width,'EdgeColor',color_1);
rectangle('Position',[pos_2(2)-0.5 pos_2(1)-0.5 1 1],'LineWidth',line_width,'EdgeColor',color_2);
rectangle('Position',[pos_3(2)-0.5 pos_3(1)-0.5 1 1],'LineWidth',line_width,'EdgeColor',color_3);

caxis([0 max(max(Elast))]);
t2=title(['E_T=',num2str(floor(E_T/(Nx*Ny)))],'FontWeight','normal');
c=colorbar;
c.Label.String='e(x,y)';
c.Ticks=[0 floor(max(max(Elast)))];
c.Label.FontSize=7;
c.FontSize=7;
c.LineWidth=line_width;
c.Color= 'k';

%COMPARISONS
%enzymes at pos_1, pos_2 and pos_3
xlist = [0.1 1 10 100];
ticks_length = 0.03;
%ax3 = subplot(2,3,4);
ax3 = subplot(3,3,4);
E_T = flux_comp(:,1);
e_pos_1 = E_prof(pos_1(1):Ny:end,pos_1(2):Ny:end);
e_pos_2 = E_prof(pos_2(1):Ny:end,pos_2(2):Ny:end);
e_pos_3 = E_prof(pos_3(1):Ny:end,pos_3(2):Ny:end);

%p_enz = plot(E_T,e_pos_1,E_T,e_pos_2,E_T,e_pos_3,'Linewidth',line_width);
p_enz = semilogx(E_T,e_pos_1,E_T,e_pos_2,E_T,e_pos_3,'Linewidth',line_width);

p_enz(1).Color = color_1;
p_enz(2).Color = color_2;
p_enz(3).Color = color_3;

set(ax3, 'YLim', [0 e_pos_max],'XLim', [0 floor(E_T_last)]);
set(gca,'Linewidth',line_width,'XAxisLocation','top','YAxisLocation','left');
xlabel('E_T');
%xticks([1 10^(floor(log10(E_T_last)))]);
xticks(xlist);
ax3.TickLength = [ticks_length,0.025];
ylabel('e_i');
%yticks([0 400 800]);
%yline(400);
%xline(200);

%inset
% %axinset = subplot(2,3,6);
% axinset = subplot(3,3,3);
% 
% inset = semilogx(E_T,e_pos_1,E_T,e_pos_2,E_T,e_pos_3,'Linewidth',line_width);
% 
% inset(1).Color = color_1;
% inset(2).Color = color_2;
% inset(3).Color = color_3;
% 
% set(axinset, 'YLim', [400 e_pos_max],'XLim', [200 floor(E_T_last)]);
% %set(gca,'Linewidth',line_width);
% %xlabel('E_T');
% %xticks([200 floor(E_T_last)]);
% xticks([]);
% %ylabel('e_i');
% %yticks([400 floor(e_pos_max)]);
% yticks([]);

%relative difference of returns at pos_1, pos_2 and pos_3
ax4 = subplot(3,3,7);
ret_pos_1 = returns(pos_1(1):Ny:end,pos_1(2):Ny:end);
ret_pos_2 = returns(pos_2(1):Ny:end,pos_2(2):Ny:end);
ret_pos_3 = returns(pos_3(1):Ny:end,pos_3(2):Ny:end);

%p_ret = plot(E_T,abs(ret_pos_2./ret_pos_1),E_T,abs(ret_pos_3./ret_pos_1),'LineWidth',line_width);
p_ret = semilogx(E_T,abs(ret_pos_2./ret_pos_1),E_T,abs(ret_pos_3./ret_pos_1),'LineWidth',line_width);

p_ret(1).Color = [0.5,0,0.5];
p_ret(2).Color = [0,0.5,0.5];

set(ax4,'YLim',[0 1], 'XLim', [0 floor(E_T_last)]);
set(gca,'Linewidth',line_width);
xlabel('E_T');
%xticks([1 10^(floor(log10(E_T_last)))]);
xticks(xlist);
ax4.TickLength = [ticks_length,0.025];
%ylabel('returns rel. change');
ylabel('returns ratio');
yticks([0 1]);
box off;

%ratio of fluxes between optimal and uniform
%ax5 = subplot(2,3,6);
ax5 = subplot(3,3,5);

J_opt = flux_comp(:,2);
J_u = flux_comp(:,3);
J_cl = flux_comp(:,4);

colororder({'k','k'})
yyaxis left
%p_fluxes_ratio = plot(E_T,J_opt./J_u,'LineWidth',line_width,'Color',[1 0.53 0]);
p_fluxes_ratio_l = semilogx(E_T,J_opt./J_u,'LineWidth',line_width);
p_fluxes_ratio_l.Color = [1 0.53 0];
ylabel('J_{opt}/J_{uniform}','Color', [1 0.53 0]);

yyaxis right
p_fluxes_ratio_r = semilogx(E_T,J_opt./J_cl,'LineWidth',line_width);
p_fluxes_ratio_r.Color = '#4a88b9';
ylabel('J_{opt}/J_{bound}', 'Color', '#4a88b9');
%ylim([1,1.5]);
%yticks([1 1.25 1.5])

%p_fluxes = plot(E_T,J_u,E_T_cl,J_cl,E_T,J_opt,'LineWidth',line_width);
%p_fluxes = loglog(E_T,J_u,E_T_cl,J_cl,E_T,J_opt,'LineWidth',line_width);
%p_fluxes(1).Color = 'k';
%p_fluxes(2).Color = '#4a88b9';
%p_fluxes(3).Color = [1 0.53 0];

set(ax5, 'XLim', [0 floor(E_T_last)]);
set(gca,'Linewidth',line_width);
xlabel('E_T');
xticks(xlist);
ax5.TickLength = [ticks_length,0.025];
box off;

% DeltaET between optimal and uniform and optimal and clustered
ax6 = subplot(3,3,8);
delta_ETu = zeros(length(J_u),1);
delta_ETcl = zeros(length(J_cl),1);
ETs = 0:0.001:floor(E_T_last);
Jopt_int = spline(E_T, J_opt, ETs)';
for i = 1:length(J_u)
    [deltaJ_min, min_i] = min(abs(J_u(i) - Jopt_int));
    delta_ETu(i) = (E_T(i)/ETs(min_i)-1)*100;
end
for i = 1:length(J_cl)
    [deltaJ_min, min_i] = min(abs(J_cl(i) - Jopt_int));
    delta_ETcl(i) = (E_T(i)/ETs(min_i)-1)*100;
end

%yyaxis left
p_deltaET_u = semilogx(E_T, delta_ETu, E_T, delta_ETcl,'LineWidth',line_width);
p_deltaET_u(1).Color = [1 0.53 0];
p_deltaET_u(2).Color = '#4a88b9';
ylabel('(\DeltaE_{T}/E_{T}) [%]');
%yticks([0 50 100 150 200])
title('extra enz. expenditure','FontWeight','normal');

%yyaxis right
%p_deltaET_cl = semilogx(E_T, delta_ETcl,'LineWidth',line_width);
%p_deltaET_cl.Color = '#4a88b9';
%ylabel('(\DeltaE_{T}/E_{T})_{bound}','Color', '#4a88b9');

set(ax6,'XLim', [0 floor(E_T_last)]);
set(gca,'Linewidth',line_width);
xlabel('E_T');
xticks(xlist);
ax6.TickLength = [ticks_length,0.025];
box off;

%FIXING THE SIZES OF THE PLOTS
%FRAMES HORIZONTAL CONFIG
%3x2 config
%ax1.Position = [0.12 0.65 0.3 0.3];
%ax1.XLabel.Position = [15 -1];
%ax2.Position = [0.55 0.65 0.3 0.3];
%ax2.XLabel.Position = [15 -1];
%c.Position = [0.87 0.65 0.02 0.3];
%c.Label.Position = [2 max(max(Elast))/2];
%ax3.Position = [0.12 0.37 0.3 0.18];
%axinset.Position = [0.2 0.44 0.1 0.08];
%ax4.Position = [0.12 0.1 0.3 0.18];
%ax5.Position = [0.55 0.37 0.3 0.18];
%ax6.Position = [0.55 0.1 0.3 0.18];

%FRAMES VERTICAL CONFIG (Sized on a figure 6.5x6.5)
%2x3 config
ax1.Position = [0.05 0.4 0.2 0.2];
ax1.XLabel.Position = [15 -1];
ax1.YLabel.Position = [-1 15];
ax2.Position = [0.05 0.1 0.2 0.2];
ax2.XLabel.Position = [15 -1];
ax2.YLabel.Position = [-1 15];
c.Position = [0.27 0.1 0.013 0.2];
c.Label.Position = [2 max(max(Elast))/2];
ax3.Position = [0.4 0.41 0.2 0.15];
% axinset.Position = [0.45 0.48 0.08 0.05];
ax4.Position = [0.4 0.13 0.2 0.15];
ax5.Position = [0.7 0.41 0.2 0.15];
ax6.Position = [0.7 0.13 0.2 0.15];

pause(1);
print('frames','-dpdf');
print('frames','-dsvg');
close gcf;
open('frames.pdf')