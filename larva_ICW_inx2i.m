clear;

Lx=1000; % length of x axis (um)
Ly=3600; % length of y axis (um)
Nx=25; % number of grid-point in x axis
Ny=90; % number of grid-point in y axis
delta_x=Lx/Nx; % size of unit grid-point in x axis (um)
delta_y=Ly/Ny; % size of unit grid-point in y axis (um)
Nt=8000; % number of time interval in the simulation

tscale=0.8; % time scale (s)
cscale=1000; % Ca ion concentration scale (uM)
K_A=0.0005; % Radio of steady concentration IP3/AKH
akhbase=240; % Basic effective AKH concentration (uM)
p_0=K_A*akhbase;% Basic IP3 concentration (uM)

sigma_c=0.02*cscale; % Fluctuation constant of Ca2+ in cytoplasm (uM)

k1 = 400; % Forward reaction rate of IP3 binding to IPR [μM^(-1) s^(-1)]
k1_minus = 52; % Backward reaction rate of IP3 binding to IPR [s^(-1)]

k2 = 20; % Forward reaction rate of Ca2+ binding to activating site of IPR [μM^(-1) s^(-1)]
k2_minus = 1.64; % Backward reaction rate of Ca2+ binding to activating site of IPR [s^(-1)]

gamma = 2; % Radio of effective volume cytoplasm/ER.
delta = 1; % Rate constant of Ca ion flux through cell membrane

Vs = 0.9*tscale*cscale; % Maximum rate of SERCA pump.[μM s^(-1)]
Ks = 0.1*cscale; % Half activation constant of Ca2+ binding to SERCA pump (μM)

Vp = 0.05*tscale*cscale; % Maximum rate of outflux through cell membrane [μM s^(-1)]
Kp = 0.3*cscale; % Half activation constant of outflux through cell membrane (μM)

alpha1 = 0.07*tscale*cscale*K_A; % IP3 contributing to the flux into the cell.[s^(-1)]
alpha0 = 0.01*tscale*cscale+akhbase*alpha1; % Constant flux into the cell.[μM s^(-1)]

kf = 1.11*tscale; % Maximal rate of Ca2+ release through IPR.[s^(-1)]
kleak = 0.02*tscale;  % Rate constant of Ca2+ leak from ER.[s^(-1)]

K1=k1_minus/k1/K_A+akhbase; % Rate constant characterizing IP3 binding to IPR.(uM)
K5=k2_minus/k2*cscale; % Rate constant characterizing Ca2+ binding to activating site of IPR.(uM)

y=0.3; % Fraction of activated IPR.

J_irp=@(c,ce,p)kf*((p+akhbase).*c.*(1-y)./(p+K1)./(c+K5)).^3.*(ce-c); % Ca2+ flux through IRP [μM s^(-1)]
J_leak=@(c,ce,p)kleak*(ce-c); % Ca2+ leak flux from ER [μM s^(-1)] 
J_serca=@(c,ce,p)Vs*c.^2./(Ks^2+c.^2); % Ca2+ flux through SERCA pump [μM s^(-1)]
J_mem=@(c,ce,p)delta*(alpha0+alpha1*p-Vp*c.^2./(Kp^2+c.^2)); % Ca2+ flux through cell membrane [μM s^(-1)] 

dc=@(c,ce,p)J_irp(c,ce,p)+J_leak(c,ce,p)-J_serca(c,ce,p)+J_mem(c,ce,p); % dynamic equation of Ca2+ in the cytoplasm
dce=@(c,ce,p)-gamma*(J_irp(c,ce,p)+J_leak(c,ce,p)-J_serca(c,ce,p)); % dynamic equation of Ca2+ in the ER

fixpoint=[0.235 ,20.6]*cscale; % Ca2+ concentration distribution in the steady state (uM)

var_c=ones(Nx,Ny,Nt)*fixpoint(1);  % Ca2+ concerntration of cytoplasm in per unit grid per time
var_ce=ones(Nx,Ny,Nt)*fixpoint(2);  % Ca2+ concerntration of ER in per unit grid per time
var_akh=zeros(Nx,Ny,Nt); % AKH concerntration in per unit grid per time

delta_t=0.1; % Time interval in the simulation (s)
D_var_c=0*tscale; % Effective diffusion constant of Ca2+ in cytoplasm [μM^2 s^(-1)]
D_var_ce=0; % Effective diffusion constant of Ca2+ in ER [μM^2 s^(-1)]
D_var_akh=150*tscale; % Effective diffusion constant of AKH [μM^2 s^(-1)]
v_flow=45*tscale; % Transport speed of AKH.[um s^(-1)]
delay=1/50*tscale; % Delay rate of AKH. [s^(-1)]
r=delta_t/delta_y*v_flow/4; 

%%% initiation
var_akh(:,1:10,1)=100;

writerObj = VideoWriter('larva_ICW_inx2i.avi'); 
open(writerObj);
x=Lx/Nx:Lx/Nx:Lx;  % for PDE not for lattice model
y=Ly/Ny:Ly/Ny:Ly;

%%% equation solving, using Crank-Nicholson method
for ti=2:1:Nt
    for xi=2:Nx-1
        %%% Store the coefficients for implicitly solving tridiagonal matrices in partial differential equations
        a=ones(1,Ny);
        b=zeros(1,Ny);
        c=r*ones(1,Ny-1);
        d=-r*ones(1,Ny-1);
        
        %%% Boundary condition in y axis
        b(1)=0;c(1)=-1;
        b(Ny)=0;d(Ny-1)=-1;
        
        for yi=2:Ny-1
            b(yi)=var_akh(xi,yi,ti-1)-r*(var_akh(xi,yi+1,ti-1)-var_akh(xi,yi-1,ti-1))...
                +delta_t*D_var_akh*((var_akh(xi-1,yi,ti-1)+var_akh(xi+1,yi,ti-1)-2*var_akh(xi,yi,ti-1))/delta_x^2+...
                (var_akh(xi,yi-1,ti-1)+var_akh(xi,yi+1,ti-1)-2*var_akh(xi,yi,ti-1))/delta_y^2)-...
                delta_t*delay*var_akh(xi,yi,ti-1);
        end

        %%% solving tridiagonal matrices equations of AKH
        var_akh(xi,:,ti)=crout(a,c,d,b);   
    end
    
    %%% Boundary condition in x axis
    var_akh(1,:,ti)=var_akh(2,:,ti);
    var_akh(Nx,:,ti)=var_akh(Nx-1,:,ti);
    
    %%% periodic AKH pulse from anterior body
    if mod(ti,3000)==0
        var_akh(:,1:10,ti)=100;
    end
    
    real_a=var_akh(:,:,ti);

    %%% Solving Ca ion dynamic equation in cytoplasm 
    for xi=2:Nx-1
        for yi=2:Ny-1
            var_c(xi,yi,ti)=var_c(xi,yi,ti-1)+delta_t*dc(var_c(xi,yi,ti-1),var_ce(xi,yi,ti-1),real_a(xi,yi))+...
                delta_t*D_var_c*((var_c(xi-1,yi,ti-1)+var_c(xi+1,yi,ti-1)-2*var_c(xi,yi,ti-1))/delta_x^2+...
                (var_c(xi,yi-1,ti-1)+var_c(xi,yi+1,ti-1)-2*var_c(xi,yi,ti-1))/delta_y^2)+sigma_c*sqrt(delta_t)*randn;
        end
    end
    
    %%% Boundary condition
    for xi=2:Nx-1
        yi=1;
        var_c(xi,yi,ti)=var_c(xi,yi,ti-1)+delta_t*dc(var_c(xi,yi,ti-1),var_ce(xi,yi,ti-1),real_a(xi,yi))+...
                delta_t*D_var_c*((var_c(xi-1,yi,ti-1)+var_c(xi+1,yi,ti-1)-2*var_c(xi,yi,ti-1))/delta_x^2+...
                2*(var_c(xi,yi+1,ti-1)-var_c(xi,yi,ti-1))/delta_y^2)+sigma_c*sqrt(delta_t)*randn;
        yi=Ny;
        var_c(xi,yi,ti)=var_c(xi,yi,ti-1)+delta_t*dc(var_c(xi,yi,ti-1),var_ce(xi,yi,ti-1),real_a(xi,yi))+...
                delta_t*D_var_c*((var_c(xi-1,yi,ti-1)+var_c(xi+1,yi,ti-1)-2*var_c(xi,yi,ti-1))/delta_x^2+...
                2*(var_c(xi,yi-1,ti-1)-var_c(xi,yi,ti-1))/delta_y^2)+sigma_c*sqrt(delta_t)*randn;
    end

    for yi=2:Ny-1
        xi=1;
        var_c(xi,yi,ti)=var_c(xi,yi,ti-1)+delta_t*dc(var_c(xi,yi,ti-1),var_ce(xi,yi,ti-1),real_a(xi,yi))+...
            delta_t*D_var_c*((var_c(xi,yi-1,ti-1)+var_c(xi,yi+1,ti-1)-2*var_c(xi,yi,ti-1))/delta_y^2+...
            2*(var_c(xi+1,yi,ti-1)-var_c(xi,yi,ti-1))/delta_x^2)+sigma_c*sqrt(delta_t)*randn;
        xi=Nx;
        var_c(xi,yi,ti)=var_c(xi,yi,ti-1)+delta_t*dc(var_c(xi,yi,ti-1),var_ce(xi,yi,ti-1),real_a(xi,yi))+...
            delta_t*D_var_c*((var_c(xi,yi-1,ti-1)+var_c(xi,yi+1,ti-1)-2*var_c(xi,yi,ti-1))/delta_y^2+...
            2*(var_c(xi-1,yi,ti-1)-var_c(xi,yi,ti-1))/delta_x^2)+sigma_c*sqrt(delta_t)*randn;
    end

    var_c(1,1,ti)=var_c(1,1,ti-1)+delta_t*dc(var_c(1,1,ti-1),var_ce(1,1,ti-1),real_a(1,1))+...
        delta_t*D_var_c*(2*(var_c(1,2,ti-1)-var_c(1,1,ti-1))/delta_y^2+...
        2*(var_c(2,1,ti-1)-var_c(1,1,ti-1))/delta_x^2)+sigma_c*sqrt(delta_t)*randn;
    var_c(1,Ny,ti)=var_c(1,Ny,ti-1)+delta_t*dc(var_c(1,Ny,ti-1),var_ce(1,Ny,ti-1),real_a(1,Ny))+...
        delta_t*D_var_c*(2*(var_c(1,Ny-1,ti-1)-var_c(1,Ny,ti-1))/delta_y^2+...
        2*(var_c(2,Ny,ti-1)-var_c(1,Ny,ti-1))/delta_x^2)+sigma_c*sqrt(delta_t)*randn;
    var_c(Nx,1,ti)=var_c(Nx,1,ti-1)+delta_t*dc(var_c(Nx,1,ti-1),var_ce(Nx,1,ti-1),real_a(Nx,1))+...
        delta_t*D_var_c*(2*(var_c(Nx,2,ti-1)-var_c(Nx,1,ti-1))/delta_y^2+...
        2*(var_c(Nx-1,1,ti-1)-var_c(Nx,1,ti-1))/delta_x^2)+sigma_c*sqrt(delta_t)*randn;
    var_c(Nx,Ny,ti)=var_c(Nx,Ny,ti-1)+delta_t*dc(var_c(Nx,Ny,ti-1),var_ce(Nx,Ny,ti-1),real_a(Nx,Ny))+...
        delta_t*D_var_c*(2*(var_c(Nx,Ny-1,ti-1)-var_c(Nx,Ny,ti-1))/delta_y^2+...
        2*(var_c(Nx-1,Ny,ti-1)-var_c(Nx,Ny,ti-1))/delta_x^2)+sigma_c*sqrt(delta_t)*randn;


    %%% Solving Ca ion dynamic equation in ER
    for xi=2:Nx-1
        for yi=2:Ny-1
            var_ce(xi,yi,ti)=var_ce(xi,yi,ti-1)+delta_t*dce(var_c(xi,yi,ti-1),var_ce(xi,yi,ti-1),real_a(xi,yi))+...
                delta_t*D_var_ce*((var_ce(xi-1,yi,ti-1)+var_ce(xi+1,yi,ti-1)-2*var_ce(xi,yi,ti-1))/delta_x^2+...
                (var_ce(xi,yi-1,ti-1)+var_ce(xi,yi+1,ti-1)-2*var_ce(xi,yi,ti-1))/delta_y^2);
        end
    end

    for xi=2:Nx-1
        yi=1;
        var_ce(xi,yi,ti)=var_ce(xi,yi,ti-1)+delta_t*dce(var_c(xi,yi,ti-1),var_ce(xi,yi,ti-1),real_a(xi,yi))+...
            delta_t*D_var_ce*((var_ce(xi-1,yi,ti-1)+var_ce(xi+1,yi,ti-1)-2*var_ce(xi,yi,ti-1))/delta_x^2+...
                2*(var_ce(xi,yi+1,ti-1)-var_ce(xi,yi,ti-1))/delta_y^2);
        yi=Ny;
        var_ce(xi,yi,ti)=var_ce(xi,yi,ti-1)+delta_t*dce(var_c(xi,yi,ti-1),var_ce(xi,yi,ti-1),real_a(xi,yi))+...
            delta_t*D_var_ce*((var_ce(xi-1,yi,ti-1)+var_ce(xi+1,yi,ti-1)-2*var_ce(xi,yi,ti-1))/delta_x^2+...
                2*(var_ce(xi,yi-1,ti-1)-var_ce(xi,yi,ti-1))/delta_y^2);
    end

    for yi=2:Ny-1
        xi=1;
        var_ce(xi,yi,ti)=var_ce(xi,yi,ti-1)+delta_t*dce(var_c(xi,yi,ti-1),var_ce(xi,yi,ti-1),real_a(xi,yi))+...
            delta_t*D_var_ce*((var_ce(xi,yi-1,ti-1)+var_ce(xi,yi+1,ti-1)-2*var_ce(xi,yi,ti-1))/delta_y^2+...
                2*(var_ce(xi+1,yi,ti-1)-var_ce(xi,yi,ti-1))/delta_x^2);
        xi=Nx;
        var_ce(xi,yi,ti)=var_ce(xi,yi,ti-1)+delta_t*dce(var_c(xi,yi,ti-1),var_ce(xi,yi,ti-1),real_a(xi,yi))+...
            delta_t*D_var_ce*((var_ce(xi,yi-1,ti-1)+var_ce(xi,yi+1,ti-1)-2*var_ce(xi,yi,ti-1))/delta_y^2+...
                2*(var_ce(xi-1,yi,ti-1)-var_ce(xi,yi,ti-1))/delta_x^2);
    end

    var_ce(1,1,ti)=var_ce(1,1,ti-1)+delta_t*dce(var_c(1,1,ti-1),var_ce(1,1,ti-1),real_a(xi,yi))+...
        delta_t*D_var_ce*(2*(var_ce(1,2,ti-1)-var_ce(1,1,ti-1))/delta_y^2+...
        2*(var_ce(2,1,ti-1)-var_ce(1,1,ti-1))/delta_x^2);
    var_ce(1,Ny,ti)=var_ce(1,Ny,ti-1)+delta_t*dce(var_c(1,Ny,ti-1),var_ce(1,Ny,ti-1),real_a(xi,yi))+...
        delta_t*D_var_ce*(2*(var_ce(1,Ny-1,ti-1)-var_ce(1,Ny,ti-1))/delta_y^2+...
        2*(var_ce(2,Ny,ti-1)-var_ce(1,Ny,ti-1))/delta_x^2);
    var_ce(Nx,1,ti)=var_ce(Nx,1,ti-1)+delta_t*dce(var_c(Nx,1,ti-1),var_ce(Nx,1,ti-1),real_a(xi,yi))+...
        delta_t*D_var_ce*(2*(var_ce(Nx,2,ti-1)-var_ce(Nx,1,ti-1))/delta_y^2+...
        2*(var_ce(Nx-1,1,ti-1)-var_ce(Nx,1,ti-1))/delta_x^2);
    var_ce(Nx,Ny,ti)=var_ce(Nx,Ny,ti-1)+delta_t*dce(var_c(Nx,Ny,ti-1),var_ce(Nx,Ny,ti-1),real_a(xi,yi))+...
        delta_t*D_var_ce*(2*(var_ce(Nx,Ny-1,ti-1)-var_ce(Nx,Ny,ti-1))/delta_y^2+...
        2*(var_ce(Nx-1,Ny,ti-1)-var_ce(Nx,Ny,ti-1))/delta_x^2);
    
    %%% draw the ICW dynamics vedio
    if mod(ti,50)==1
        imagesc(x, y, var_c(:,:,ti)'); 
        daspect([1 1 1]);
        load Newcolormap.mat;
        colormap(NewColormap); 
        colorbar;
        clim([0,2.5]*cscale);
        set(gca,'xtick',[],'xticklabel',[])
        set(gca,'ytick',[],'yticklabel',[])
        title(['t=',num2str((ti-1)*delta_t),'s']);
        writeVideo(writerObj, getframe(gcf));
    end
    ti
end

close(writerObj);

%% Frames of simulated global ICWs in the fat body
load Newcolormap.mat;
colormap(NewColormap);
t_frame=[3010,3310,3610,3910,4210,4510];
for i=1:6
    subplot(1,6,i);
    ti=t_frame(i);
    imagesc(x, y, var_c(:,:,ti)'); 
    % imagesc(x, y, var_akh(:,:,ti)'); 
    daspect([1 1 1]);
    set(gca,'xtick',[],'xticklabel',[])
    set(gca,'ytick',[],'yticklabel',[])
    clim([0,2.5]*cscale);
    axis off;
end

%% Frames of Akh distribution dynamics
load Newcolormap.mat;
colormap(NewColormap);
t_frame=[3010,3310,3610,3910,4210,4510];
for i=1:6
    subplot(1,6,i);
    ti=t_frame(i);
    imagesc(x, y, var_akh(:,:,ti)'); 
    daspect([1 1 1]);
    set(gca,'xtick',[],'xticklabel',[])
    set(gca,'ytick',[],'yticklabel',[])
    clim([0,100]);
    axis off;
end

%% Colormap of the tempo-spacial ICW
t=delta_t*(1:Nt);
yt_track=reshape(var_c(round(Nx/2),:,:),Ny,Nt);

load Newcolormap.mat;
colormap(NewColormap);
imagesc(t,y,yt_track);
clim([0,2.5]*cscale);
set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])

%% Colormap of the tempo-spacial Akh dynamics
t=delta_t*(1:Nt);
yt_track=reshape(var_akh(round(Nx/2),:,:),Ny,Nt);

load Newcolormap.mat;
colormap(NewColormap);
imagesc(t,y,yt_track);
clim([0,100]);
set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])

%% Calculation of ICW propagation speed
t=delta_t*(1:Nt);
yt_track=reshape(var_c(round(Nx/2),:,:),Ny,Nt);

%%% record of the displacemnet of ICW in y axis at every moment
peak_start=[];
ps_check=[];
peak_end=[];
pe_check=[];
for i=2:Nt-1
    if yt_track(30,i)>yt_track(30,i-1) && yt_track(30,i)>yt_track(30,i+1) && yt_track(30,i)>400
        peak_start=[peak_start,i*delta_t];
        ps_check=[ps_check,yt_track(30,i)];
    end
    if yt_track(Ny-30,i)>yt_track(Ny-30,i-1) && yt_track(Ny-30,i)>yt_track(Ny-30,i+1) && yt_track(Ny-30,i)>400
        peak_end=[peak_end,i*delta_t];
        pe_check=[pe_check,yt_track(Ny-30,i)];
    end
end

% ICW speed (um/s)
wave_num=length(peak_end);
wave_speed=(Ny-30-30)*delta_y./(peak_end(1:wave_num)-peak_start(1:wave_num));
mean_speed=mean(wave_speed(2:end));

%% Calculation of half-peak wave width
Ny_eff=Ny-9-30;
yt_track_single=reshape(var_c(round(Nx/2),30:Ny-10,3001:6000),Ny_eff,3000);

t_halfpeak=zeros(1,Ny_eff);
for yj=1:Ny_eff
    t_start=0;
    t_end=0;    
    try
        t_start=find(yt_track_single(yj,:)>400,1,"first");
        t_end=find(yt_track_single(yj,:)>400,1,'last');
        t_halfpeak(yj)=(t_end-t_start)*delta_t;
    end
end
t_halfpeak(t_halfpeak==0)=[];

mean_t_hf=mean(t_halfpeak);
mean_pk=mean(max(yt_track_single));

%% Analysis of coorelation along the ICW propagation
Ny_cor=floor((Ny-1)/3); % every 3 rows as a analytical unit
cor_rela=zeros(Ny_cor,Nt); % Record the square of the variance/mean at each time, at each y coordinate
mean_rela=zeros(Ny_cor,Nt); % Record the mean at each time, at each y coordinate
for ti=1:Nt
    for yi=1:Ny_cor
        var_tem=var(var_c(:,3*yi-1:3*yi+1,ti),0,"all");
        mean_tem=mean(var_c(:,3*yi-1:3*yi+1,ti),"all");
        cor_rela(yi,ti)=var_tem/mean_tem^2;
        mean_rela(yi,ti)=mean_tem;
    end
end

t=delta_t*(1:Nt); 
y1=3*delta_y*(2:Ny_cor);

%%% Colormap which shows tempo-spacial dynamics of the square of CV coefficient 
figure;
imagesc(t,y1,cor_rela(2:end,:));
colorbar;
clim([0,0.5]);
xlabel('time/s')
ylabel('y coordinate/um')

%%% Plot of the square of CV coefficient at the wavefront
tstart_cor=70;
tend_cor=120;
kt_y=(tend_cor-tstart_cor)/Ny_cor;
mid_cor=zeros(Ny_cor-1,1);
tplt_cor=zeros(Ny_cor-1,1);
for yi=2:Ny_cor
    ti=round((tstart_cor+yi*kt_y)/delta_t);
    tplt_cor(yi-1)=tstart_cor+yi*kt_y;
    mid_cor(yi-1)=cor_rela(yi,ti);
end

figure;
plot(tplt_cor,mid_cor,'LineWidth',2);
xlabel('time/s')
ylabel('square of CV coefficient')
ylim([0,0.5]);