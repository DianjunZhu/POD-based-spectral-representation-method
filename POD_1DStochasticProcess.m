clc;clear;
num=1000; % Number of stochastic process
T=102.4;
detaT=0.1;
t=0:detaT:T-detaT; % Time
nt=round(T/detaT);
dw=2*pi/T;
fw=pi/detaT;
w=0:dw:fw-dw; % Frequency
%% EPSD
[tt,ww]=meshgrid(t,w);
A0=10;
b=0.075;
A=A0-b.*tt;
Sn=A.^3./4.*ww.^2.*exp(-A.*abs(ww));
vars=sum(Sn)*dw*2; % Standard deviation
%% Generation of stochastic process
ncut=6;% The truncation number of POD
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
ylim([0, 5]);
zlim([0, 1.5]);
xticks([0:25:T])
yticks([0:1:5])
zticks([0:0.5:1.5])
view([51.4,31.2]);
title('Target EPSD')
set(gca,'FontName','Times New Roman','FontSize',12);
grid on
%
figure
subplot(2,1,1);
hold on
plot(t,h(1,:),'-','LineWidth',1.0)
plot(t,h(2,:),'-','LineWidth',1.0)
plot(t,h(3,:),'-','LineWidth',1.0)
plot(t,h(4,:),'-','LineWidth',1.0)
plot(t,h(5,:),'-','LineWidth',1.0)
plot(t,h(6,:),'-','LineWidth',1.0)
xlabel('Time (s)','FontSize',12)
ylabel('Principal coordinate','FontSize',12)
legend({'Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6'},'Location','best')
set(gca,'FontName','Times New Roman','FontSize',12);
legend('FontName','Times New Roman','FontSize',12);
xlim([0, T]);
ylim([-4,4]);
xticks([0:20:T])
yticks([-4:2:4])
grid on
box on
subplot(2,1,2);
hold on
plot(w,a(:,1),'-','LineWidth',1.0)
plot(w,a(:,2),'-','LineWidth',1.0)
plot(w,a(:,3),'-','LineWidth',1.0)
plot(w,a(:,4),'-','LineWidth',1.0)
plot(w,a(:,5),'-','LineWidth',1.0)
plot(w,a(:,6),'-','LineWidth',1.0)
xlabel('Frequency (rad/s)','FontSize',12)
ylabel('Eigenvector','FontSize',12)
set(gca,'FontName','Times New Roman','FontSize',12);
xlim([0, 5]);
ylim([-0.4,0.4]);
xticks([0:1:5])
yticks([-0.4:0.2:0.4])
grid on
box on
%
figure
hold on
p1=plot(w,Sn(:,1/detaT)/vars(1/detaT),'-b','LineWidth',2);
p2=plot(w,Sn(:,50/detaT)/vars(50/detaT),'-b','LineWidth',2);
p3=plot(w,Sn(:,100/detaT)/vars(100/detaT),'-b','LineWidth',2);
p4=plot(w,Gs(:,1/detaT).^2/vars(1/detaT),'--r','LineWidth',2);
p5=plot(w,Gs(:,50/detaT).^2/vars(50/detaT),'--r','LineWidth',2);
p6=plot(w,Gs(:,100/detaT).^2/vars(100/detaT),'--r','LineWidth',2);
xlabel('Frequency (rad/s)','FontSize',12)
ylabel('Standardized EPSD','FontSize',12)
legend([p1 p4],{'Target EPSD','POD'},'Location','northeast')
set(gca,'FontName','Times New Roman','FontSize',12);
legend('FontName','Times New Roman','FontSize',12);
xlim([0, 3]);
ylim([0, 1.5]);
xticks([0:1:3])
yticks([0:0.5:1.5])
grid on
box on
%
index=randi([1, num], 2, 1);
figure
subplot(2,1,2);
acc_plot1=plot(t,Um(:,index(1)),'-b', 'LineWidth', 1.5 );
xlim([0, T]);
ylim([-3,3]);
ylabel('Acceleration (cm/s^2)','FontSize',12)
xlabel('Time (s)','FontSize',12)
xticks([0:20:T])
yticks([-3:2:3])
title(['Sample ',num2str(index(1))])
set(gca,'FontName','Times New Roman','FontSize',12);
grid on
box on
subplot(2,1,1);
hold on
acc_plot2=plot(t,Um(:,index(2)),'-b', 'LineWidth', 1.5 );
xlim([0, T]);
ylim([-3,3]);
ylabel('Acceleration (cm/s^2)','FontSize',12)
xlabel('Time (s)','FontSize',12)
xticks([0:20:T])
yticks([-3:2:3])
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
ylim([-1,1]);
xticks([0:20:T])
yticks([-1:0.5:1])
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
ylim([0,2]);
xticks([0:20:T])
yticks([0:0.5:2])
grid on
box on