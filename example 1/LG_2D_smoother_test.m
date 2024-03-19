clear all; clc; clf;
currentFolder = pwd;
addpath(genpath(currentFolder));
% Model parameters;
m.ss.Talfa=1; % Sampling time interval;
m.ss.qw=0.5; % System noise spectual density;
m.ss.C=[0.75,0.75];
m.ss.A=[1,m.ss.Talfa;0,1]; % State transition matrix;
m.ss.Q = m.ss.qw.*[m.ss.Talfa^3/3,m.ss.Talfa^2/2;m.ss.Talfa^2/2,m.ss.Talfa];
m.ss.Qnn=m.ss.Q(1,1); m.ss.Qnl=m.ss.Q(1,2); m.ss.gamaln=m.ss.Qnl'/m.ss.Qnn; m.ss.Qll=m.ss.Q(2,2);m.ss.Qllbar=m.ss.Qll-m.ss.gamaln*m.ss.Qnn*m.ss.gamaln';
m.ss.gamanl=m.ss.gamaln';
m.ss.M=m.ss.Q*m.ss.C'; m.ss.gamany=m.ss.M(1,1)'/m.ss.Qnn;
m.ss.gamalybar=(m.ss.M(2,1)-m.ss.gamaln*m.ss.M(1,1))/m.ss.Qllbar;
m.ss.An=m.ss.A(1,2); m.ss.Al=m.ss.A(2,2);m.ss.Albar=m.ss.Al-m.ss.gamaln*m.ss.An;
m.ss.H=[1,0];
m.ss.Hbar1=-m.ss.gamany'*m.ss.An-m.ss.gamalybar'*m.ss.Albar; m.ss.Hbar2=m.ss.H(1,2)+m.ss.gamalybar';
m.ss.Hbar=[m.ss.Hbar1,m.ss.Hbar2];
m.ss.R0=1;
m.ss.R = m.ss.R0 + m.ss.C*m.ss.Q*m.ss.C';


m.ss.D1=m.ss.Hbar(1,1)+m.ss.Hbar(1,2)*m.ss.Albar;
m.ss.D2=m.ss.Hbar(1,2)-m.ss.gamalybar';
m.ss.D3=m.ss.D2*m.ss.gamaln+m.ss.gamany';

m.ss.m0 = [0;1]; m.ss.mn0=m.ss.m0(1,1); m.ss.ml0=m.ss.m0(2,1);
m.ss.P0 = eye(2); m.ss.Pn0=m.ss.P0(1,1); m.ss.Pl0=m.ss.P0(2,2);
m.ss.dimx=2; m.ss.dimy=1; m.ss.dimxn=1; m.ss.dimxl=1;

m.ss.N=100; % Number of filtering particles
m.ss.Nthres=5*m.ss.N/5;
m.ss.Ns=100;
m.ss.rep=10;
m.ss.T = 101; % Length of data
K=1000; % Number of MC runs
GKFdep=cell(1,K);
GBPFdep=cell(1,K);
GBRBPFind=cell(1,K);
GBRBPFdep=cell(1,K);

GPSdep=cell(1,K);
GRBPSind=cell(1,K);
GRBPSdep=cell(1,K);
GMHRBPSdep=cell(1,K);
GMHRBPPdep=cell(1,K);
GKSdep=cell(1,K);

XnkfdeperrorArr=zeros(K,m.ss.T);     XlkfdeperrorArr=zeros(K,m.ss.T);
XnbpfdeperrorArr=zeros(K,m.ss.T);     XlbpfdeperrorArr=zeros(K,m.ss.T);
XnbrbpfinderrorArr=zeros(K,m.ss.T);   XlbrbpfinderrorArr=zeros(K,m.ss.T);
XnbrbpfdeperrorArr=zeros(K,m.ss.T);   XlbrbpfdeperrorArr=zeros(K,m.ss.T);

XnpsdeperrorArr=zeros(K,m.ss.T);      XlpsdeperrorArr=zeros(K,m.ss.T);
XnrbpsinderrorArr=zeros(K,m.ss.T);    XlrbpsinderrorArr=zeros(K,m.ss.T);
XnrbpsdeperrorArr=zeros(K,m.ss.T);    XlrbpsdeperrorArr=zeros(K,m.ss.T);
XnmhrbpsdeperrorArr=zeros(K,m.ss.T);    XlmhrbpsdeperrorArr=zeros(K,m.ss.T);
XnmhrbppdeperrorArr=zeros(K,m.ss.T);    XlmhrbppdeperrorArr=zeros(K,m.ss.T);
XnksdeperrorArr=zeros(K,m.ss.T);       XlksdeperrorArr=zeros(K,m.ss.T);


Timekfdep=0;Timebpfdep=0;Timebrbpfind=0;Timebrbpfdep=0;
Timepsdep=0;Timerbpsind=0;Timerbpsdep=0;Timemhrbpsdep=0;Timemhrbppdep=0;Timeksdep=0;
% i=1;
% while i<=K
%     while 1
%         try
%             ind=1;
for i=1:K
    disp(['Current Monte Carlo iteration : ', num2str(i)]);
    % Simulate the system;
    Xtrue=zeros(m.ss.dimx,m.ss.T);
    y=zeros(m.ss.dimy,m.ss.T);
    for k=1:m.ss.T
        if k==1
            Xtrue(:,k)=mvnrnd(m.ss.m0,m.ss.P0,1)';
        else
            w=sqrtm(m.ss.Q)*randn(m.ss.dimx,1);
            Xtrue(:,k)=m.ss.A*Xtrue(:,k-1)+w;
            v=m.ss.C*w+sqrtm(m.ss.R0)*randn(1);
            y(k)=m.ss.H*Xtrue(:,k)+v;%k时刻观测向量；
        end
    end
    m.x=Xtrue;
    m.y=y;     % Store the measurements
    %% Filters;
    GKFdep{i}=KF_LG_2D_asyn_cor(m);
    GBPFdep{i}=BPF_LG_2D_asyn_cor(m);
    GBRBPFind{i}=BRBPF_LG_2D_ind(m);
    GBRBPFdep{i}=BRBPF_LG_2D_asyn_cor(m);
    %% Smoothers;
    GPSdep{i}= ffbsi_asyn_cor(GBPFdep{i},m);
    GRBPSind{i}=RBPS_ind(GBRBPFind{i},m);
    GRBPSdep{i}=RBPS_LG2D_asyn_cor(GBRBPFdep{i},m);
    GMHRBPSdep{i}=MH_RBPS_LG2D_asyn_cor(GBRBPFdep{i},m);
    GMHRBPPdep{i}=MH_RBPSnew_LG2D_asyn_cor(GBRBPFdep{i},m);
    GKSdep{i}=KS_asyn_cor(GKFdep{i},m);
    %         catch
    %             ind=0;
    %         end
    %         if ind==1
    %             break;
    %         end
    %     end
    %% Filtering errors;
    XnkfdeperrorArr(i,:)=GKFdep{i}.Xerror(1,:);
    XlkfdeperrorArr(i,:)=GKFdep{i}.Xerror(2,:);
    XnbpfdeperrorArr(i,:)=GBPFdep{i}.Xerror(1,:);
    XlbpfdeperrorArr(i,:)=GBPFdep{i}.Xerror(2,:);
    XnbrbpfinderrorArr(i,:)=GBRBPFind{i}.Xerror(1,:);
    XlbrbpfinderrorArr(i,:)=GBRBPFind{i}.Xerror(2,:);
    XnbrbpfdeperrorArr(i,:)=GBRBPFdep{i}.Xerror(1,:);
    XlbrbpfdeperrorArr(i,:)=GBRBPFdep{i}.Xerror(2,:);
    %% Smoothing errors;
    XnpsdeperrorArr(i,:)=GPSdep{i}.Xerror(1,:);
    XlpsdeperrorArr(i,:)=GPSdep{i}.Xerror(2,:);
    XnrbpsinderrorArr(i,:)=GRBPSind{i}.Xerror(1,:);
    XlrbpsinderrorArr(i,:)=GRBPSind{i}.Xerror(2,:);
    XnrbpsdeperrorArr(i,:)=GRBPSdep{i}.Xerror(1,:);
    XlrbpsdeperrorArr(i,:)=GRBPSdep{i}.Xerror(2,:);
    XnmhrbpsdeperrorArr(i,:)=GMHRBPSdep{i}.Xerror(1,:);
    XlmhrbpsdeperrorArr(i,:)=GMHRBPSdep{i}.Xerror(2,:);
    XnmhrbppdeperrorArr(i,:)=GMHRBPPdep{i}.Xerror(1,:);
    XlmhrbppdeperrorArr(i,:)=GMHRBPPdep{i}.Xerror(2,:);
    XnksdeperrorArr(i,:)=GKSdep{i}.Xkserror(1,:);
    XlksdeperrorArr(i,:)=GKSdep{i}.Xkserror(2,:);
    %% Filtering running time;
    Timekfdep=Timekfdep+GKFdep{i}.time;
    Timebpfdep=Timebpfdep+GBPFdep{i}.time;
    Timebrbpfind=Timebrbpfind+GBRBPFind{i}.time;
    Timebrbpfdep=Timebrbpfdep+GBRBPFdep{i}.time;
    %% Smoothing running time;
    Timepsdep=Timepsdep+GPSdep{i}.time;
    Timerbpsind=Timerbpsind+GRBPSind{i}.time;
    Timerbpsdep=Timerbpsdep+GRBPSdep{i}.time;
    Timemhrbpsdep=Timemhrbpsdep+GMHRBPSdep{i}.time;
    Timemhrbppdep=Timemhrbppdep+GMHRBPPdep{i}.time;
    Timeksdep=Timeksdep+GKSdep{i}.time;
%     disp(['Successful Monte Carlo iteration: ', num2str(i)]);
%     i=i+1;
end
KF_dep_RMSE_xp=sqrt(mean(XnkfdeperrorArr.^2));KF_dep_RMSE_xk=sqrt(mean(XlkfdeperrorArr.^2));
BPF_dep_RMSE_xp=sqrt(mean(XnbpfdeperrorArr.^2));BPF_dep_RMSE_xk=sqrt(mean(XlbpfdeperrorArr.^2));
BRBPF_ind_RMSE_xp=sqrt(mean(XnbrbpfinderrorArr.^2));BRBPF_ind_RMSE_xk=sqrt(mean(XlbrbpfinderrorArr.^2));
BRBPF_dep_RMSE_xp=sqrt(mean(XnbrbpfdeperrorArr.^2));BRBPF_dep_RMSE_xk=sqrt(mean(XlbrbpfdeperrorArr.^2));

PS_dep_RMSE_xp=sqrt(mean(XnpsdeperrorArr.^2));PS_dep_RMSE_xk=sqrt(mean(XlpsdeperrorArr.^2));
RBPS_ind_RMSE_xp=sqrt(mean(XnrbpsinderrorArr.^2));RBPS_ind_RMSE_xk=sqrt(mean(XlrbpsinderrorArr.^2));
RBPS_dep_RMSE_xp=sqrt(mean(XnrbpsdeperrorArr.^2));RBPS_dep_RMSE_xk=sqrt(mean(XlrbpsdeperrorArr.^2));
MHRBPS_dep_RMSE_xp=sqrt(mean(XnmhrbpsdeperrorArr.^2));MHRBPS_dep_RMSE_xk=sqrt(mean(XlmhrbpsdeperrorArr.^2));
MHRBPP_dep_RMSE_xp=sqrt(mean(XnmhrbppdeperrorArr.^2));MHRBPP_dep_RMSE_xk=sqrt(mean(XlmhrbppdeperrorArr.^2));
KS_dep_RMSE_xp=sqrt(mean(XnksdeperrorArr.^2));KS_dep_RMSE_xk=sqrt(mean(XlksdeperrorArr.^2));

Position_RMSE_Arr=[BRBPF_dep_RMSE_xp;PS_dep_RMSE_xp;RBPS_ind_RMSE_xp;RBPS_dep_RMSE_xp;MHRBPS_dep_RMSE_xp;MHRBPP_dep_RMSE_xp;KS_dep_RMSE_xp];
Velocity_RMSE_Arr=[BRBPF_dep_RMSE_xk;PS_dep_RMSE_xk;RBPS_ind_RMSE_xk;RBPS_dep_RMSE_xk;MHRBPS_dep_RMSE_xk;MHRBPP_dep_RMSE_xk;KS_dep_RMSE_xk];
Averaged_time=[Timebrbpfdep,Timepsdep,Timerbpsind,Timerbpsdep,Timemhrbpsdep,Timemhrbppdep,Timeksdep]./K;
t=0:m.ss.T-1;
figure(1);
plot(t,BRBPF_dep_RMSE_xp,'-sr',t,PS_dep_RMSE_xp,'-og',t,RBPS_ind_RMSE_xp,'-hm',t,RBPS_dep_RMSE_xp,'-<b',t,MHRBPS_dep_RMSE_xp,'-cd',t,MHRBPP_dep_RMSE_xp,'-+k',t,KS_dep_RMSE_xp,'-py','linewidth',1);
legend('D-RBPF','D-PS','I-RBPS','D-RBPS','D-MH-RBBS','D-MH-RBBP','D-KS');
xlabel('Time step');
ylabel('Position RMSE');
figure(2);
plot(t,BRBPF_dep_RMSE_xk,'-sr',t,PS_dep_RMSE_xk,'-og',t,RBPS_ind_RMSE_xk,'-hm',t,RBPS_dep_RMSE_xk,'-<b',t,MHRBPS_dep_RMSE_xk,'-cd',t,MHRBPP_dep_RMSE_xk,'-+k',t,KS_dep_RMSE_xk,'-py','linewidth',1);
legend('D-RBPF','D-PS','I-RBPS','D-RBPS','D-MH-RBBS','D-MH-RBBP','D-KS');
xlabel('Time step');
ylabel('Velocity RMSE');
disp('Time averaged position RMSEs for DRBPF, DPS, IRBPS, DRBPS, DMHRBPS, DMHRBBP, DKS are:')
mean(Position_RMSE_Arr,2)'
disp('Time averaged velocity RMSEs for DRBPF, DPS, IRBPS, DRBPS, DMHRBPS, DMHRBBP, DKS are:')
mean(Velocity_RMSE_Arr,2)'
disp('The averaged running time for  DRBPF, DPS, IRBPS, DRBPS, DMHRBPS, DMHRBBP, DKS are:')
Averaged_time





