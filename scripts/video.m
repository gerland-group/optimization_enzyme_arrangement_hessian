load '../E_prof.out';
load '../returns.out';
load '../flux_comp.out';
clf;

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
writerObj = VideoWriter('optimization','Uncompressed AVI');
writerObj.FrameRate=30;
open(writerObj);
Elast=E_prof((m-(Ny-1)):m,1:n);

figure('pos',[0 0 1000 400],'Resize','off');
colormap('bone');
colormap(flipud(colormap));

set(0,'DefaultAxesFontSize',12)

set(gca,'FontName','Helvetica','FontSize',12,'nextplot','replace');
set(gcf,'Renderer','opengl');
set(gcf,'color','w');
%MACROS FOR STYLE
line_width=1.5;
x_com_size=0.4;
y_com_size=0.4;
line_style_sources='--';
%positions to compare [row column]
pos_1=[10 9];
color_1 = 'b';
pos_2=[25 24];
color_2 = 'r';

e_pos_1_last = Elast(pos_1(1),pos_1(2));
e_pos_2_last = Elast(pos_2(1),pos_2(2));
e_pos_max = max(e_pos_1_last,  e_pos_2_last);
E_T_last = sum(Elast(:))/(Nx*Ny);
J_comp_uniform_max = max((flux_comp(:,2)./flux_comp(:,3)));
J_comp_uniform_min = min((flux_comp(:,2)./flux_comp(:,3)));
%J_comp_cluster_max = max((flux_comp(:,2)./flux_comp(:,4)));
%J_comp_cluster_min = min((flux_comp(:,2)./flux_comp(:,4)));


e_pos_1=[];
e_pos_2=[];
ret_pos_1=[];
ret_pos_2=[];
J_u=[];
J_opt=[];
J_cl=[];
E_T=[];
%video
k=10*Ny;
for i=1:k:m
    %SYSTEM FRAME
    ax1 = subplot(2,3,[1 4]);
    E_frame=E_prof(i:i+(Ny-1),1:n);
    R_frame=returns(i:i+(Ny-1),1:n);
    
    imagesc(E_frame);
    xlabel('x');
    xticks([0.5 Nx+0.5]);
    xticklabels({0,Nx});
    set (gca,'Ydir','normal','Linewidth',line_width);
    ylabel('y');
    yticks([0.5 Ny+0.5]);
    yticklabels({0,Ny});

    E_T_i = sum(E_frame(:))/(Nx*Ny);
    
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

    if E_T_i>0
        [y_grid,x_grid] = ndgrid(1:size(E_frame,1),1:size(E_frame,2));
        y_com = sum(y_grid(:).*E_frame(:))/(E_T_i*Nx*Ny);
        x_com = sum(x_grid(:).*E_frame(:))/(E_T_i*Nx*Ny);
        
        rectangle('Position',[x_com-x_com_size/2 y_com-y_com_size/2 x_com_size y_com_size],'Curvature',[1 1],'Linewidth',line_width, 'FaceColor','w');
    end
    
    rectangle('Position',[0.5 0.5 Nx Ny],'Linewidth',line_width);
    rectangle('Position',[pos_1(2)-0.5 pos_1(1)-0.5 1 1],'LineWidth',line_width,'EdgeColor',color_1);
    rectangle('Position',[pos_2(2)-0.5 pos_2(1)-0.5 1 1],'LineWidth',line_width,'EdgeColor',color_2);
    
    caxis([0 max(max(Elast))]);
    c=colorbar;
    c.Label.String='e(x,y)';
    c.Ticks=[0 floor(max(max(Elast)))];
    c.Label.FontSize=12;
    c.FontSize=12;
    c.LineWidth=line_width;
    c.Color= 'k';
    %COMPARISONS
    %enzymes at pos_1 and pos_2
    ax2 = subplot(2,3,2);
    E_T = [E_T E_T_i];
    e_pos_1_i = E_frame(pos_1(1),pos_1(2));
    e_pos_1 = [e_pos_1,e_pos_1_i];
    e_pos_2_i = E_frame(pos_2(1),pos_2(2));
    e_pos_2 = [e_pos_2,e_pos_2_i];
    
    p_enz = plot(E_T,e_pos_1,E_T,e_pos_2,'Linewidth',line_width);
    p_enz(1).Color = color_1;
    p_enz(2).Color = color_2;
    
    set(ax2, 'YLim', [0 e_pos_max],'XLim', [0 floor(E_T_last)]);
    set(gca,'Linewidth',line_width);
    xlabel('E_T');
    xticks([0 floor(E_T_last)]);
    ylabel('e_1, e_2');
    yticks([0 floor(e_pos_max)]);
    
    %returns at pos_1 and pos_2
    ax3 = subplot(2,3,5);
    ret_pos_1_i = R_frame(pos_1(1),pos_1(2));
    ret_pos_1 = [ret_pos_1,ret_pos_1_i];
    ret_pos_2_i = R_frame(pos_2(1),pos_2(2));
    ret_pos_2 = [ret_pos_2,ret_pos_2_i];
    
    %p_ret = plot(E_T,ret_pos_1,E_T,ret_pos_2,'Linewidth',line_width);
    %p_ret(1).Color = color_1;
    %p_ret(2).Color = color_2;
    
    p_ret = plot(E_T,abs(1-ret_pos_2./ret_pos_1),'LineWidth',line_width,'Color',[0.5,0,0.5]);
    
    set(ax3,'YLim',[0 1], 'XLim', [0 floor(E_T_last)]);
    set(gca,'Linewidth',line_width);
    xlabel('E_T');
    xticks([0 floor(E_T_last)]);
    ylabel('returns rel. change');
    yticks([0 1]);
    
    %relative difference of fluxes between optimal and uniform
    ax4 = subplot(2,3,3);
    J_opt_i = flux_comp((i-1)/k+1,2);
    J_opt = [J_opt J_opt_i];
    J_u_i = flux_comp((i-1)/k+1,3);
    J_u = [J_u J_u_i];
    
    p_fluxes = plot(E_T,J_opt./J_u,'LineWidth',line_width,'Color',[0 0.5 0]);
    
    set(ax4,'YLim',[0 floor(J_comp_uniform_max+1)], 'XLim', [0 floor(E_T_last)]);
    set(gca,'Linewidth',line_width);
    xlabel('E_T');
    xticks([0 floor(E_T_last)]);
    ylabel('J_{opt}/J_u');
    %yticks([0 4]);
    %YTickLabels = num2str(log10([10^2,10^3,10^4]), '10^%d');
    
    %relative difference of fluxes between optimal and cluster
%     ax5 = subplot(2,3,6);
%     J_cl_i = flux_comp((i-1)/k+1,4);
%     J_cl = [J_cl J_cl_i];
%     
%     p_fluxes = plot(E_T,J_opt./J_cl,'LineWidth',line_width,'Color',[0 0.5 0]);
%     
%     set(ax5,'YLim',[0 floor(J_comp_cluster_max+1)], 'XLim', [0 floor(E_T_last)]);
%     set(gca,'Linewidth',line_width);
%     xlabel('E_T');
%     xticks([0 floor(E_T_last)]);
%     ylabel('J_{opt}/J_{cl}');
    %yticks([0 1]);
    %YTickLabels = num2str(log10([10^2,10^3,10^4]), '10^%d');
    
    %REGULATING SUBPLOTS POSITIONS AND RATIOS
    ax1.Position = [0.05 0.15 0.32 0.8];
    ax2.Position = [0.55 0.65 0.15 0.3];
    ax3.Position = [0.55 0.15 0.15 0.3];
    ax4.Position = [0.8 0.65 0.15 0.3];
%     ax5.Position = [0.8 0.15 0.15 0.3];
    
    %GETFRAME
    frame=getframe(gcf);
    writeVideo(writerObj,frame);
        
end
close(writerObj);

