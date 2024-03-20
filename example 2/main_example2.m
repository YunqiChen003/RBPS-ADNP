clear all; clc; clf;
currentFolder = pwd;
addpath(genpath(currentFolder));
m.ss.delta=1;
m.ss.o=-pi/60;
m.ss.m0=[2000,2000,10,10]';m.ss.mn0=m.ss.m0(1:2);m.ss.ml0=m.ss.m0(3:end);
m.ss.P0=diag([4,4,4,4]);m.ss.Pn0=m.ss.P0(1:2,1:2);m.ss.Pl0=m.ss.P0(3:end,3:end);

m.ss.Qxx=0.25*[m.ss.delta^3/3*eye(2),m.ss.delta^2/2*eye(2);m.ss.delta^2/2*eye(2),eye(2)];
% m.ss.Qxx=[2*(m.ss.o*m.ss.delta-sin(m.ss.o*m.ss.delta))/m.ss.o^3,0,(1-cos(m.ss.o*m.ss.delta))/m.ss.o^2,(m.ss.o*m.ss.delta-sin(m.ss.o*m.ss.delta))/m.ss.o^2;...
%      0,2*(m.ss.o*m.ss.delta-sin(m.ss.o*m.ss.delta))/m.ss.o^3,-(m.ss.o*m.ss.delta-sin(m.ss.o*m.ss.delta))/m.ss.o^2,(1-cos(m.ss.o*m.ss.delta))/m.ss.o^2;...
%      (1-cos(m.ss.o*m.ss.delta))/m.ss.o^2,-(m.ss.o*m.ss.delta-sin(m.ss.o*m.ss.delta))/m.ss.o^2,m.ss.delta,0;...
%      (m.ss.o*m.ss.delta-sin(m.ss.o*m.ss.delta))/m.ss.o^2,(1-cos(m.ss.o*m.ss.delta))/m.ss.o^2,0,m.ss.delta];
m.ss.Qnn=m.ss.Qxx(1:2,1:2);m.ss.Qnl=m.ss.Qxx(1:2,3:end);m.ss.Qll=m.ss.Qxx(3:end,3:end);
m.ss.R=diag([300,10^-6]);
m.ss.Qxy=[3,3,5,5;0,0,0,0]';m.ss.Qny=m.ss.Qxy(1:2,:);m.ss.Qly=m.ss.Qxy(3:end,:);
m.ss.Sigma=[m.ss.Qxx,m.ss.Qxy;m.ss.Qxy',m.ss.R];

m.ss.A=[1,0,sin(m.ss.o*m.ss.delta)/m.ss.o,(cos(m.ss.o*m.ss.delta)-1)/m.ss.o;0,1,(1-cos(m.ss.o*m.ss.delta))/m.ss.o,sin(m.ss.o*m.ss.delta)/m.ss.o;...
    0,0,cos(m.ss.o*m.ss.delta),-sin(m.ss.o*m.ss.delta);0,0,sin(m.ss.o*m.ss.delta),cos(m.ss.o*m.ss.delta)]; m.ss.An=m.ss.A(1:2,3:end);
m.ss.gamanl=m.ss.Qnn\m.ss.Qnl;
m.ss.Al=m.ss.A(3:end,3:end);
m.ss.Albar=m.ss.Al-m.ss.gamanl'*m.ss.An;
m.ss.Qllbar=m.ss.Qll-m.ss.gamanl'*m.ss.Qnl;
m.ss.Qlybar=m.ss.Qly-m.ss.gamanl'*m.ss.Qny;

m.ss.gamany=m.ss.Qnn\m.ss.Qny; m.ss.gamalybar=m.ss.Qllbar\m.ss.Qlybar;
m.ss.R0=m.ss.R-m.ss.gamany'*m.ss.Qnn*m.ss.gamany-m.ss.gamalybar'*m.ss.Qllbar*m.ss.gamalybar;
m.ss.Hbar1=-m.ss.gamany'*m.ss.An-m.ss.gamalybar'*m.ss.Albar; m.ss.Hbar2=m.ss.gamalybar';
m.ss.H=zeros(2,2);
m.ss.Hbar=[m.ss.Hbar1,m.ss.Hbar2];
m.ss.D1=m.ss.Hbar1+m.ss.Hbar2*m.ss.Albar; m.ss.D2=m.ss.Hbar2-m.ss.gamalybar'; m.ss.D3=m.ss.D2*m.ss.gamanl'+m.ss.gamany';

m.ss.N=500;m.ss.Ns=100;
m.ss.dimx=4;m.ss.dimxn=2;m.ss.dimxl=2;
m.ss.dimy=2;
m.ss.Nthres=5*m.ss.N/5;
m.ss.rep=30;
K=500;m.ss.T=101;

GBPFdep=cell(1,K);
GBRBPFind=cell(1,K);
GBRBPFdep=cell(1,K);

GPSdep=cell(1,K);
GRBPSind=cell(1,K);
GRBPSdep=cell(1,K);
GMHRBPSdep=cell(1,K);
GMHRBPPdep=cell(1,K);


XbrbpfdeperrorArr=zeros(m.ss.dimx,m.ss.T,K);
XpsdeperrorArr=zeros(m.ss.dimx,m.ss.T,K);
XrbpsinderrorArr=zeros(m.ss.dimx,m.ss.T,K);
XrbpsdepmferrorArr=zeros(m.ss.dimx,m.ss.T,K);
XmhrbpsdepmferrorArr=zeros(m.ss.dimx,m.ss.T,K);
XmhrbppdepmferrorArr=zeros(m.ss.dimx,m.ss.T,K);
XrbpsdeprtserrorArr=zeros(m.ss.dimx,m.ss.T,K);
XmhrbpsdeprtserrorArr=zeros(m.ss.dimx,m.ss.T,K);
XmhrbppdeprtserrorArr=zeros(m.ss.dimx,m.ss.T,K);


Timebrbpfdep=0;
Timepsdep=0;Timerbpsind=0;
Timerbpsdepmf=0;Timemhrbpsdepmf=0;Timemhrbppdepmf=0;
Timerbpsdeprts=0;Timemhrbpsdeprts=0;Timemhrbppdeprts=0;
 
for i=1:K
    disp(['Current Monte Carlo iteration : ', num2str(i)]);
    Xtrue=zeros(m.ss.dimx,m.ss.T);
    y=zeros(m.ss.dimy,m.ss.T);
    for k=1:m.ss.T
        if k==1
            Xtrue(:,k)=mvnrnd(m.ss.m0,m.ss.P0,1)';
        else
            noise=(mvnrnd(zeros(m.ss.dimx+m.ss.dimy,1),m.ss.Sigma,1))';
            Xtrue(:,k)=sys_eq(Xtrue(:,k-1),m)+noise(1:m.ss.dimx,:);
            y(:,k)=mea_eq(Xtrue(:,k))+noise(m.ss.dimx+1:end,:);%k时刻观测向量；
        end
    end
    m.x=Xtrue;
    m.y=y;     % Store the measurements
    %% Filters;
    GBPFdep{i}=BPF_asyn_cor(m);
    GBRBPFind{i}=SORBPF_ind(m);
    GBRBPFdep{i}=SORBPF_asyn_cor(m);
    %% Smoothers;
    GPSdep{i}= ffbsi_asyn_cor(GBPFdep{i},m);
    GRBPSind{i}=RBPS_ind(GBRBPFind{i},m);
    GRBPSdep{i}=RBPS_asyn_cor(GBRBPFdep{i},m);
    GMHRBPSdep{i}=MH_RBPS_asyn_cor(GBRBPFdep{i},m);
    GMHRBPPdep{i}=MH_RBPSnew_asyn_cor(GBRBPFdep{i},m);
    %% Filtering errors;
    XbrbpfdeperrorArr(:,:,i)=GBRBPFdep{i}.Xerror;
    %% Smoothing errors;
    XpsdeperrorArr(:,:,i)=GPSdep{i}.Xerror;
    XrbpsinderrorArr(:,:,i)=GRBPSind{i}.Xerror;
    XrbpsdepmferrorArr(:,:,i)=GRBPSdep{i}.Xerror1;
    XrbpsdeprtserrorArr(:,:,i)=GRBPSdep{i}.Xerror2;
    XmhrbpsdepmferrorArr(:,:,i)=GMHRBPSdep{i}.Xerror1;
    XmhrbpsdeprtserrorArr(:,:,i)=GMHRBPSdep{i}.Xerror2;
    XmhrbppdepmferrorArr(:,:,i)=GMHRBPPdep{i}.Xerror1;
    XmhrbppdeprtserrorArr(:,:,i)=GMHRBPPdep{i}.Xerror2;
    %% Filtering running time;
    Timebrbpfdep=Timebrbpfdep+GBRBPFdep{i}.time;
     %% Smoothing running time;
    Timepsdep=Timepsdep+GPSdep{i}.time;
    Timerbpsind=Timerbpsind+GRBPSind{i}.time;
    Timerbpsdepmf=Timerbpsdepmf+GRBPSdep{i}.time1;
    Timemhrbpsdepmf=Timemhrbpsdepmf+GMHRBPSdep{i}.time1;
    Timemhrbppdepmf=Timemhrbppdepmf+GMHRBPPdep{i}.time1;
    Timerbpsdeprts=Timerbpsdeprts+GRBPSdep{i}.time2;
    Timemhrbpsdeprts=Timemhrbpsdeprts+GMHRBPSdep{i}.time2;
    Timemhrbppdeprts=Timemhrbppdeprts+GMHRBPPdep{i}.time2;
end

BRBPFdep_RMSEposition=sqrt(mean(XbrbpfdeperrorArr(1,:,:).^2+XbrbpfdeperrorArr(2,:,:).^2,3));
BRBPFdep_RMSEvelocity=sqrt(mean(XbrbpfdeperrorArr(3,:,:).^2+XbrbpfdeperrorArr(4,:,:).^2,3));
BRBPFdep_RMSE=sqrt(mean(XbrbpfdeperrorArr.^2,3));

PS_dep_RMSEposition=sqrt(mean(XpsdeperrorArr(1,:,:).^2+XpsdeperrorArr(2,:,:).^2,3));
PS_dep_RMSEvelocity=sqrt(mean(XpsdeperrorArr(3,:,:).^2+XpsdeperrorArr(4,:,:).^2,3));
PS_dep_RMSE=sqrt(mean(XpsdeperrorArr.^2,3));


RBPS_ind_RMSEposition=sqrt(mean(XrbpsinderrorArr(1,:,:).^2+XrbpsinderrorArr(2,:,:).^2,3));
RBPS_ind_RMSEvelocity=sqrt(mean(XrbpsinderrorArr(3,:,:).^2+XrbpsinderrorArr(4,:,:).^2,3));
RBPS_ind_RMSE=sqrt(mean(XrbpsinderrorArr.^2,3));

RBPS_dep_mf_RMSEposition=sqrt(mean(XrbpsdepmferrorArr(1,:,:).^2+XrbpsdepmferrorArr(2,:,:).^2,3));
RBPS_dep_mf_RMSEvelocity=sqrt(mean(XrbpsdepmferrorArr(3,:,:).^2+XrbpsdepmferrorArr(4,:,:).^2,3));
RBPS_dep_mf_RMSE=sqrt(mean(XrbpsdepmferrorArr.^2,3));

RBPS_dep_rts_RMSEposition=sqrt(mean(XrbpsdeprtserrorArr(1,:,:).^2+XrbpsdeprtserrorArr(2,:,:).^2,3));
RBPS_dep_rts_RMSEvelocity=sqrt(mean(XrbpsdeprtserrorArr(3,:,:).^2+XrbpsdeprtserrorArr(4,:,:).^2,3));
RBPS_dep_rts_RMSE=sqrt(mean(XrbpsdeprtserrorArr.^2,3));

MHRBPS_dep_mf_RMSEposition=sqrt(mean(XmhrbpsdepmferrorArr(1,:,:).^2+XmhrbpsdepmferrorArr(2,:,:).^2,3));
MHRBPS_dep_mf_RMSEvelocity=sqrt(mean(XmhrbpsdepmferrorArr(3,:,:).^2+XmhrbpsdepmferrorArr(4,:,:).^2,3));
MHRBPS_dep_mf_RMSE=sqrt(mean(XmhrbpsdepmferrorArr.^2,3));

MHRBPS_dep_rts_RMSEposition=sqrt(mean(XmhrbpsdeprtserrorArr(1,:,:).^2+XmhrbpsdeprtserrorArr(2,:,:).^2,3));
MHRBPS_dep_rts_RMSEvelocity=sqrt(mean(XmhrbpsdeprtserrorArr(3,:,:).^2+XmhrbpsdeprtserrorArr(4,:,:).^2,3));
MHRBPS_dep_rts_RMSE=sqrt(mean(XmhrbpsdeprtserrorArr.^2,3));

MHRBPP_dep_mf_RMSEposition=sqrt(mean(XmhrbppdepmferrorArr(1,:,:).^2+XmhrbppdepmferrorArr(2,:,:).^2,3));
MHRBPP_dep_mf_RMSEvelocity=sqrt(mean(XmhrbppdepmferrorArr(3,:,:).^2+XmhrbppdepmferrorArr(4,:,:).^2,3));
MHRBPP_dep_mf_RMSE=sqrt(mean(XmhrbppdepmferrorArr.^2,3));

MHRBPP_dep_rts_RMSEposition=sqrt(mean(XmhrbppdeprtserrorArr(1,:,:).^2+XmhrbppdeprtserrorArr(2,:,:).^2,3));
MHRBPP_dep_rts_RMSEvelocity=sqrt(mean(XmhrbppdeprtserrorArr(3,:,:).^2+XmhrbppdeprtserrorArr(4,:,:).^2,3));
MHRBPP_dep_rts_RMSE=sqrt(mean(XmhrbppdeprtserrorArr.^2,3));

Position_RMSE_Arr=[BRBPFdep_RMSEposition;PS_dep_RMSEposition;RBPS_ind_RMSEposition;...
                   RBPS_dep_mf_RMSEposition;MHRBPS_dep_mf_RMSEposition;MHRBPP_dep_mf_RMSEposition;...
                   RBPS_dep_rts_RMSEposition;MHRBPS_dep_rts_RMSEposition;MHRBPP_dep_rts_RMSEposition];
Velocity_RMSE_Arr=[BRBPFdep_RMSEvelocity;PS_dep_RMSEvelocity;RBPS_ind_RMSEvelocity;...
                   RBPS_dep_mf_RMSEvelocity;MHRBPS_dep_mf_RMSEvelocity;MHRBPP_dep_mf_RMSEvelocity;...
                   RBPS_dep_rts_RMSEvelocity;MHRBPS_dep_rts_RMSEvelocity;MHRBPP_dep_rts_RMSEvelocity];
Px_RMSE_Arr=[BRBPFdep_RMSE(1,:);PS_dep_RMSE(1,:);RBPS_ind_RMSE(1,:);...
             RBPS_dep_mf_RMSE(1,:);MHRBPS_dep_mf_RMSE(1,:);MHRBPP_dep_mf_RMSE(1,:);...
             RBPS_dep_rts_RMSE(1,:);MHRBPS_dep_rts_RMSE(1,:);MHRBPP_dep_rts_RMSE(1,:)];
Py_RMSE_Arr=[BRBPFdep_RMSE(2,:);PS_dep_RMSE(2,:);RBPS_ind_RMSE(2,:);...
             RBPS_dep_mf_RMSE(2,:);MHRBPS_dep_mf_RMSE(2,:);MHRBPP_dep_mf_RMSE(2,:);...
             RBPS_dep_rts_RMSE(2,:);MHRBPS_dep_rts_RMSE(2,:);MHRBPP_dep_rts_RMSE(2,:)];
Vx_RMSE_Arr=[BRBPFdep_RMSE(3,:);PS_dep_RMSE(3,:);RBPS_ind_RMSE(3,:);...
             RBPS_dep_mf_RMSE(3,:);MHRBPS_dep_mf_RMSE(3,:);MHRBPP_dep_mf_RMSE(3,:);...
             RBPS_dep_rts_RMSE(3,:);MHRBPS_dep_rts_RMSE(3,:);MHRBPP_dep_rts_RMSE(3,:)];
Vy_RMSE_Arr=[BRBPFdep_RMSE(4,:);PS_dep_RMSE(4,:);RBPS_ind_RMSE(4,:);...
             RBPS_dep_mf_RMSE(4,:);MHRBPS_dep_mf_RMSE(4,:);MHRBPP_dep_mf_RMSE(4,:);...
             RBPS_dep_rts_RMSE(4,:);MHRBPS_dep_rts_RMSE(4,:);MHRBPP_dep_rts_RMSE(4,:)];
Averaged_time=[Timebrbpfdep,Timepsdep,Timerbpsind,...
               Timerbpsdepmf,Timemhrbpsdepmf,Timemhrbppdepmf,...
               Timerbpsdeprts,Timemhrbpsdeprts,Timemhrbppdeprts]./K;

t=0:m.ss.T-1;
figure(1);
subplot(2,1,1)
plot(t,BRBPFdep_RMSEposition,'-r',t,PS_dep_RMSEposition,'-g',t,RBPS_ind_RMSEposition,'-m',t,RBPS_dep_mf_RMSEposition,'-b',t,MHRBPS_dep_mf_RMSEposition,'-c',t,MHRBPP_dep_mf_RMSEposition,'-k','linewidth',1);
legend('D-RBPF','D-PS','I-RBPS','D-RBPS','D-MH-RBPS1','D-MH-RBPS2');
xlabel('Time (s)');
ylabel('Position RMSE (m)');
subplot(2,1,2)
plot(t,RBPS_dep_mf_RMSEposition,'-b',t,MHRBPS_dep_mf_RMSEposition,'-c',t,MHRBPP_dep_mf_RMSEposition,'-k',t,RBPS_dep_rts_RMSEposition,'-r',t,MHRBPS_dep_rts_RMSEposition,'-g',t,MHRBPP_dep_rts_RMSEposition,'-m','linewidth',1);
legend('D-RBPS','D-MH-RBPS1','D-MH-RBPS2','D-RBPS*','D-MH-RBPS1*','D-MH-RBPS2*');
xlabel('Time (s)');
ylabel('Position RMSE (m)');
figure(2);
subplot(2,1,1)
plot(t,BRBPFdep_RMSEvelocity,'-r',t,PS_dep_RMSEvelocity,'-g',t,RBPS_ind_RMSEvelocity,'-m',t,RBPS_dep_mf_RMSEvelocity,'-b',t,MHRBPS_dep_mf_RMSEvelocity,'-c',t,MHRBPP_dep_mf_RMSEvelocity,'-k','linewidth',1);
legend('D-RBPF','D-PS','I-RBPS','D-RBPS','D-MH-RBPS1','D-MH-RBPS2');
xlabel('Time (s)');
ylabel('Velocity RMSE (m/s)');
subplot(2,1,2)
plot(t,RBPS_dep_mf_RMSEvelocity,'-b',t,MHRBPS_dep_mf_RMSEvelocity,'-c',t,MHRBPP_dep_mf_RMSEvelocity,'-k',t,RBPS_dep_rts_RMSEvelocity,'-r',t,MHRBPS_dep_rts_RMSEvelocity,'-g',t,MHRBPP_dep_rts_RMSEvelocity,'-m','linewidth',1);
legend('D-RBPS','D-MH-RBPS1','D-MH-RBPS2','D-RBPS*','D-MH-RBPS1*','D-MH-RBPS2*');
xlabel('Time (s)');
ylabel('Velocity RMSE (m)');
figure(3);
subplot(2,1,1)
plot(t,BRBPFdep_RMSE(1,:),'-r',t,PS_dep_RMSE(1,:),'-g',t,RBPS_ind_RMSE(1,:),'-m',t,RBPS_dep_mf_RMSE(1,:),'-b',t,MHRBPS_dep_mf_RMSE(1,:),'-c',t,MHRBPP_dep_mf_RMSE(1,:),'-k','linewidth',1);
legend('D-RBPF','D-PS','I-RBPS','D-RBPS','D-MH-RBPS1','D-MH-RBPS2');
xlabel('Time (s)');
ylabel('$p_{k}^{x}$ RMSE (m)','Interpreter', 'latex');
subplot(2,1,2)
plot(t,RBPS_dep_mf_RMSE(1,:),'-b',t,MHRBPS_dep_mf_RMSE(1,:),'-c',t,MHRBPP_dep_mf_RMSE(1,:),'-k',t,RBPS_dep_rts_RMSE(1,:),'-r',t,MHRBPS_dep_rts_RMSE(1,:),'-g',t,MHRBPP_dep_rts_RMSE(1,:),'-m','linewidth',1);
legend('D-RBPS','D-MH-RBPS1','D-MH-RBPS2','D-RBPS*','D-MH-RBPS1*','D-MH-RBPS2*');
xlabel('Time (s)');
ylabel('$p_{k}^{x}$ RMSE (m)','Interpreter', 'latex');
figure(4);
subplot(2,1,1)
plot(t,BRBPFdep_RMSE(2,:),'-r',t,PS_dep_RMSE(2,:),'-g',t,RBPS_ind_RMSE(2,:),'-m',t,RBPS_dep_mf_RMSE(2,:),'-b',t,MHRBPS_dep_mf_RMSE(2,:),'-c',t,MHRBPP_dep_mf_RMSE(2,:),'-k','linewidth',1);
legend('D-RBPF','D-PS','I-RBPS','D-RBPS','D-MH-RBPS1','D-MH-RBPS2');
xlabel('Time (s)');
ylabel('$p_{k}^{y}$ RMSE (m)','Interpreter', 'latex');
subplot(2,1,2)
plot(t,RBPS_dep_mf_RMSE(2,:),'-b',t,MHRBPS_dep_mf_RMSE(2,:),'-c',t,MHRBPP_dep_mf_RMSE(2,:),'-k',t,RBPS_dep_rts_RMSE(2,:),'-r',t,MHRBPS_dep_rts_RMSE(2,:),'-g',t,MHRBPP_dep_rts_RMSE(2,:),'-m','linewidth',1);
legend('D-RBPS','D-MH-RBPS1','D-MH-RBPS2','D-RBPS*','D-MH-RBPS1*','D-MH-RBPS2*');
xlabel('Time (s)');
ylabel('$p_{k}^{y}$ RMSE (m)','Interpreter', 'latex');
figure(5);
subplot(2,1,1)
plot(t,BRBPFdep_RMSE(3,:),'-r',t,PS_dep_RMSE(3,:),'-g',t,RBPS_ind_RMSE(3,:),'-m',t,RBPS_dep_mf_RMSE(3,:),'-b',t,MHRBPS_dep_mf_RMSE(3,:),'-c',t,MHRBPP_dep_mf_RMSE(3,:),'-k','linewidth',1);
legend('D-RBPF','D-PS','I-RBPS','D-RBPS','D-MH-RBPS1','D-MH-RBPS2');
xlabel('Time (s)');
ylabel('$\dot{p}_{k}^{x}$ RMSE (m/s)','Interpreter', 'latex');
subplot(2,1,2)
plot(t,RBPS_dep_mf_RMSE(3,:),'-b',t,MHRBPS_dep_mf_RMSE(3,:),'-c',t,MHRBPP_dep_mf_RMSE(3,:),'-k',t,RBPS_dep_rts_RMSE(3,:),'-r',t,MHRBPS_dep_rts_RMSE(3,:),'-g',t,MHRBPP_dep_rts_RMSE(3,:),'-m','linewidth',1);
legend('D-RBPS','D-MH-RBPS1','D-MH-RBPS2','D-RBPS*','D-MH-RBPS1*','D-MH-RBPS2*');
xlabel('Time (s)');
ylabel('$\dot{p}_{k}^{x}$ RMSE (m/s)','Interpreter', 'latex');
figure(6);
subplot(2,1,1)
plot(t,BRBPFdep_RMSE(4,:),'-r',t,PS_dep_RMSE(4,:),'-g',t,RBPS_ind_RMSE(4,:),'-m',t,RBPS_dep_mf_RMSE(4,:),'-b',t,MHRBPS_dep_mf_RMSE(4,:),'-c',t,MHRBPP_dep_mf_RMSE(4,:),'-k','linewidth',1);
legend('D-RBPF','D-PS','I-RBPS','D-RBPS','D-MH-RBPS1','D-MH-RBPS2');
xlabel('Time (s)');
ylabel('$\dot{p}_{k}^{y}$ RMSE (m/s)','Interpreter', 'latex');
subplot(2,1,2)
plot(t,RBPS_dep_mf_RMSE(4,:),'-b',t,MHRBPS_dep_mf_RMSE(4,:),'-c',t,MHRBPP_dep_mf_RMSE(4,:),'-k',t,RBPS_dep_rts_RMSE(4,:),'-r',t,MHRBPS_dep_rts_RMSE(4,:),'-g',t,MHRBPP_dep_rts_RMSE(4,:),'-m','linewidth',1);
legend('D-RBPS','D-MH-RBPS1','D-MH-RBPS2','D-RBPS*','D-MH-RBPS1*','D-MH-RBPS2*');
xlabel('Time (s)');
ylabel('$\dot{p}_{k}^{y}$ RMSE (m/s)','Interpreter', 'latex');



disp('Pos_ARMSEs for D-RBPF,D-PS,I-RBPS,D-RBPS,D-MHRBPS1,D-MHRBPS2,D-RBPS*,D-MHRBPS1*,D-MHRBPS2* are:')
mean(Position_RMSE_Arr,2)'
disp('Vel_ARMSEs for D-RBPF,D-PS,I-RBPS,D-RBPS,D-MHRBPS1,D-MHRBPS2,D-RBPS*,D-MHRBPS1*,D-MHRBPS2* are:')
mean(Velocity_RMSE_Arr,2)'
disp('Px_ARMSEs for D-RBPF,D-PS,I-RBPS,D-RBPS,D-MHRBPS1,D-MHRBPS2,D-RBPS*,D-MHRBPS1*,D-MHRBPS2* are:')
mean(Px_RMSE_Arr,2)'
disp('Py_ARMSEs for D-RBPF,D-PS,I-RBPS,D-RBPS,D-MHRBPS1,D-MHRBPS2,D-RBPS*,D-MHRBPS1*,D-MHRBPS2* are:')
mean(Py_RMSE_Arr,2)'
disp('Vx_ARMSEs for D-RBPF,D-PS,I-RBPS,D-RBPS,D-MHRBPS1,D-MHRBPS2,D-RBPS*,D-MHRBPS1*,D-MHRBPS2* are:')
mean(Vx_RMSE_Arr,2)'
disp('Vy_ARMSEs for D-RBPF,D-PS,I-RBPS,D-RBPS,D-MHRBPS1,D-MHRBPS2,D-RBPS*,D-MHRBPS1*,D-MHRBPS2* are:')
mean(Vy_RMSE_Arr,2)'
disp('The averaged running time for D-RBPF,D-PS,I-RBPS,D-RBPS,D-MHRBPS1,D-MHRBPS2,D-RBPS*,D-MHRBPS1*,D-MHRBPS2* are:')
Averaged_time





