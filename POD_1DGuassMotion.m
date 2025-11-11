clc;clear;
num=1000; % Number of stochastic ground motions
T=20.48;
detaT=0.01;
t=detaT:detaT:T; % Time
nt=round(T/detaT);
dw=2*pi/T;
fw=pi/detaT;
w=0:dw:fw-dw; % Frequency
%% Clough-Penzien power spectrum model parameters
w0=22.5;
yeta0=0.8;
gama=2.35;
nw=round(fw/dw);
[tt,ww]=meshgrid(t,w);
a=12;
b=0.65;
wg=w0-a*t./T;
yetag=yeta0+b*t./T;
wf=0.1*wg;
yetaf=yetag;
amax=196;
%% EPSD
% Modulation functions
c=6; d=3;
%1
A=[tt/c.*exp(1-tt/c)].^d;
%2
% tt1=0.8;tt2=7;Lamb=0.35;
% A(:,1:tt1/detaT)=(tt(:,1:tt1/detaT)/tt1).^2;
% A(:,tt1/detaT+1:tt2/detaT)=1;
% A(:,tt2/detaT+1:T/detaT)=exp(-Lamb.*(tt(:,tt2/detaT+1:T/detaT)-tt2));
%3
% nn0=0.15;
% A=exp(-nn0.*ww.*tt/wg./T).*[tt/c.*exp(1-tt/c)].^d;
% Generalized Clough-Penzien power spectrum model
S0=amax.^2./(gama.^2*pi.*wg.*(2*yetag+1./(2*yetag)));
S=(wg.^4+4*(yetag.^2).*(wg.^2).*(ww.^2))./(((ww.^2)-(wg.^2)).^2+4*(yetag.^2).*(wg.^2).*(ww.^2))...
    .*ww.^4./(((ww.^2)-(wf.^2)).^2+4*(yetaf.^2).*(wf.^2).*(ww.^2)).*S0;
Sn=(A.*A).*S;
vars=sum(Sn)*dw*2; % Standard deviation
%% Generation of non-stationary stochastic ground motions
ncut=3;% The truncation number of POD
randn('state',1);rand('state',1);
[ Um, a, ~, h ] = POD_fft( Sn, ncut, dw, num);% POD-FFT
% Relative error
Gs=a*h;
err=norm(sqrt(Sn)-Gs,"fro")/norm(sqrt(Sn),"fro");
disp(['The relative error of POD is: ',num2str(err)]);
%% Figures
figure
surf_EPSD=surf(t,w,Sn);
surf_EPSD.FaceAlpha=0.6;
surf_EPSD.LineStyle='-';
surf_EPSD.EdgeColor = 'flat';
colormap jet
grid on
xlabel('Time (s)','FontSize',12)
ylabel('Frequency (rad/s)','FontSize',12)
zlabel('EPSD','FontSize',12)
xlim([0, T]);
ylim([0, 120]);
zlim([0, 80]);
xticks([0:5:T])
yticks([0:20:120])
zticks([0:20:80])
view([51.4,31.2]);
title('Target EPSD')
set(gca,'FontName','Times New Roman','FontSize',12);
grid on
%
figure
subplot(2,1,1);
hold on
plot(t,h(1,:),'b-','LineWidth',1.0)
plot(t,h(2,:),'g-','LineWidth',1.0)
plot(t,h(3,:),'k-','LineWidth',1.0)
xlabel('Time (s)','FontSize',12)
ylabel('Principal coordinate','FontSize',12)
legend({'Mode 1','Mode 2','Mode 3'},'Location','northeast')
set(gca,'FontName','Times New Roman','FontSize',12);
legend('FontName','Times New Roman','FontSize',12);
xlim([0, T]);
ylim([-25,100]);
xticks([0:5:T])
yticks([-25:25:125])
grid on
box on
subplot(2,1,2);
hold on
plot(w,a(:,1),'b-','LineWidth',1.0)
plot(w,a(:,2),'g-','LineWidth',1.0)
plot(w,a(:,3),'k-','LineWidth',1.0)
xlabel('Frequency (rad/s)','FontSize',12)
ylabel('Eigenvector','FontSize',12)
set(gca,'FontName','Times New Roman','FontSize',12);
xlim([0, 120]);
ylim([-0.2,0.2]);
xticks([0:20:120])
yticks([-0.2:0.1:0.2])
grid on
box on
%
figure
hold on
p1=plot(w,100*Sn(:,1/detaT)/vars(1/detaT),'-b','LineWidth',2);
p2=plot(w,100*Sn(:,10/detaT)/vars(10/detaT),'-b','LineWidth',2);
p3=plot(w,100*Sn(:,20/detaT)/vars(20/detaT),'-b','LineWidth',2);
p4=plot(w,100*Gs(:,1/detaT).^2/vars(1/detaT),'--r','LineWidth',2);
p5=plot(w,100*Gs(:,10/detaT).^2/vars(10/detaT),'--r','LineWidth',2);
p6=plot(w,100*Gs(:,20/detaT).^2/vars(20/detaT),'--r','LineWidth',2);
xlabel('Frequency (rad/s)','FontSize',12)
ylabel('Standardized EPSD (10^-2)','FontSize',12)
legend([p1 p4],{'Target EPSD','POD'},'Location','northeast')
set(gca,'FontName','Times New Roman','FontSize',12);
legend('FontName','Times New Roman','FontSize',12);
xlim([0, 120]);
ylim([0, 1.2]);
xticks([0:20:120])
yticks([0:0.3:1.2])
grid on
box on
%
index=randi([1, num], 2, 1);
figure
subplot(2,1,2);
acc_plot1=plot(t,Um(:,index(1)),'-b', 'LineWidth', 1.5 );
xlim([0, T]);
ylim([-200,200]);
ylabel('Acceleration (cm/s^2)','FontSize',12)
xlabel('Time (s)','FontSize',12)
xticks([0:5:T])
yticks([-200:100:200])
title(['Sample ',num2str(index(1))])
set(gca,'FontName','Times New Roman','FontSize',12);
grid on
box on
subplot(2,1,1);
hold on
acc_plot2=plot(t,Um(:,index(2)),'-b', 'LineWidth', 1.5 );
xlim([0, T]);
ylim([-200,200]);
ylabel('Acceleration (cm/s^2)','FontSize',12)
xlabel('Time (s)','FontSize',12)
xticks([0:5:T])
yticks([-200:100:200])
title(['Sample ',num2str(index(2))])
set(gca,'FontName','Times New Roman','FontSize',12);
grid on
box on
%
figure
subplot(2,1,1);
hold on
plot(t,zeros(nt,1),'k-','LineWidth',1.5)
plot(t,mean(Um'),'r-','LineWidth',1.0)
xlabel('Time (s)','FontSize',12)
ylabel('Mean','FontSize',12)
legend({'Target value','MCS'},'Location','northeast')
set(gca,'FontName','Times New Roman','FontSize',12);
legend('FontName','Times New Roman','FontSize',12);
xlim([0, T]);
ylim([-20,20]);
xticks([0:5:T])
yticks([-20:10:20])
grid on
box on
subplot(2,1,2);
hold on
plot(t,sqrt(vars),'k-','LineWidth',1.5)
plot(t,std(Um'),'r-','LineWidth',1.0)
xlabel('Time (s)','FontSize',12)
ylabel('S.D.','FontSize',12)
legend({'Target value','MCS'},'Location','northeast')
set(gca,'FontName','Times New Roman','FontSize',12);
legend('FontName','Times New Roman','FontSize',12);
xlim([0, T]);
ylim([0,80]);
xticks([0:5:T])
yticks([0:20:80])
grid on
box on