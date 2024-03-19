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
K=100; % Number of MC runs
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

XnkfdeperrorArr=zeros(K,m.ss.T);          XlkfdeperrorArr=zeros(K,m.ss.T);
XnbpfdeperrorArr=zeros(K,m.ss.T);         XlbpfdeperrorArr=zeros(K,m.ss.T);
XnbrbpfinderrorArr=zeros(K,m.ss.T);       XlbrbpfinderrorArr=zeros(K,m.ss.T);
XnbrbpfdeperrorArr=zeros(K,m.ss.T);       XlbrbpfdeperrorArr=zeros(K,m.ss.T);

XnpsdeperrorArr=zeros(K,m.ss.T);          XlpsdeperrorArr=zeros(K,m.ss.T);
XnrbpsinderrorArr=zeros(K,m.ss.T);        XlrbpsinderrorArr=zeros(K,m.ss.T);

XnrbpsdepmferrorArr=zeros(K,m.ss.T);      XlrbpsdepmferrorArr=zeros(K,m.ss.T);
XnrbpsdeprtserrorArr=zeros(K,m.ss.T);     XlrbpsdeprtserrorArr=zeros(K,m.ss.T);

XnmhrbpsdepmferrorArr=zeros(K,m.ss.T);    XlmhrbpsdepmferrorArr=zeros(K,m.ss.T);
XnmhrbpsdeprtserrorArr=zeros(K,m.ss.T);   XlmhrbpsdeprtserrorArr=zeros(K,m.ss.T);

XnmhrbppdepmferrorArr=zeros(K,m.ss.T);    XlmhrbppdepmferrorArr=zeros(K,m.ss.T);
XnmhrbppdeprtserrorArr=zeros(K,m.ss.T);   XlmhrbppdeprtserrorArr=zeros(K,m.ss.T);

XnksdeperrorArr=zeros(K,m.ss.T);          XlksdeperrorArr=zeros(K,m.ss.T);


Timekfdep=0;Timebpfdep=0;Timebrbpfind=0;Timebrbpfdep=0;
Timepsdep=0;Timerbpsind=0;
Timerbpsdepmf=0;Timemhrbpsdepmf=0;Timemhrbppdepmf=0;
Timerbpsdeprts=0;Timemhrbpsdeprts=0;Timemhrbppdeprts=0;
Timeksdep=0;
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
    XnrbpsdepmferrorArr(i,:)=GRBPSdep{i}.Xerror1(1,:);
    XlrbpsdepmferrorArr(i,:)=GRBPSdep{i}.Xerror1(2,:);
    XnrbpsdeprtserrorArr(i,:)=GRBPSdep{i}.Xerror2(1,:);
    XlrbpsdeprtserrorArr(i,:)=GRBPSdep{i}.Xerror2(2,:); 
    XnmhrbpsdepmferrorArr(i,:)=GMHRBPSdep{i}.Xerror1(1,:);
    XlmhrbpsdepmferrorArr(i,:)=GMHRBPSdep{i}.Xerror1(2,:);
    XnmhrbpsdeprtserrorArr(i,:)=GMHRBPSdep{i}.Xerror2(1,:);
    XlmhrbpsdeprtserrorArr(i,:)=GMHRBPSdep{i}.Xerror2(2,:); 
    XnmhrbppdepmferrorArr(i,:)=GMHRBPPdep{i}.Xerror1(1,:);
    XlmhrbppdepmferrorArr(i,:)=GMHRBPPdep{i}.Xerror1(2,:);
    XnmhrbppdeprtserrorArr(i,:)=GMHRBPPdep{i}.Xerror2(1,:);
    XlmhrbppdeprtserrorArr(i,:)=GMHRBPPdep{i}.Xerror2(2,:);
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
    Timerbpsdepmf=Timerbpsdepmf+GRBPSdep{i}.time1;
    Timemhrbpsdepmf=Timemhrbpsdepmf+GMHRBPSdep{i}.time1;
    Timemhrbppdepmf=Timemhrbppdepmf+GMHRBPPdep{i}.time1;
    Timerbpsdeprts=Timerbpsdeprts+GRBPSdep{i}.time2;
    Timemhrbpsdeprts=Timemhrbpsdeprts+GMHRBPSdep{i}.time2;
    Timemhrbppdeprts=Timemhrbppdeprts+GMHRBPPdep{i}.time2;
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
RBPS_dep_mf_RMSE_xp=sqrt(mean(XnrbpsdepmferrorArr.^2));RBPS_dep_mf_RMSE_xk=sqrt(mean(XlrbpsdepmferrorArr.^2));
RBPS_dep_rts_RMSE_xp=sqrt(mean(XnrbpsdeprtserrorArr.^2));RBPS_dep_rts_RMSE_xk=sqrt(mean(XlrbpsdeprtserrorArr.^2));
MHRBPS_dep_mf_RMSE_xp=sqrt(mean(XnmhrbpsdepmferrorArr.^2));MHRBPS_dep_mf_RMSE_xk=sqrt(mean(XlmhrbpsdepmferrorArr.^2));
MHRBPS_dep_rts_RMSE_xp=sqrt(mean(XnmhrbpsdeprtserrorArr.^2));MHRBPS_dep_rts_RMSE_xk=sqrt(mean(XlmhrbpsdeprtserrorArr.^2));
MHRBPP_dep_mf_RMSE_xp=sqrt(mean(XnmhrbppdepmferrorArr.^2));MHRBPP_dep_mf_RMSE_xk=sqrt(mean(XlmhrbppdepmferrorArr.^2));
MHRBPP_dep_rts_RMSE_xp=sqrt(mean(XnmhrbppdeprtserrorArr.^2));MHRBPP_dep_rts_RMSE_xk=sqrt(mean(XlmhrbppdeprtserrorArr.^2));
KS_dep_RMSE_xp=sqrt(mean(XnksdeperrorArr.^2));KS_dep_RMSE_xk=sqrt(mean(XlksdeperrorArr.^2));

Position_RMSE_Arr=[BRBPF_dep_RMSE_xp;PS_dep_RMSE_xp;RBPS_ind_RMSE_xp;...
                   RBPS_dep_mf_RMSE_xp;MHRBPS_dep_mf_RMSE_xp;MHRBPP_dep_mf_RMSE_xp;...
                   RBPS_dep_rts_RMSE_xp;MHRBPS_dep_rts_RMSE_xp;MHRBPP_dep_rts_RMSE_xp;...
                   KS_dep_RMSE_xp];
Velocity_RMSE_Arr=[BRBPF_dep_RMSE_xk;PS_dep_RMSE_xk;RBPS_ind_RMSE_xk;...
                   RBPS_dep_mf_RMSE_xk;MHRBPS_dep_mf_RMSE_xk;MHRBPP_dep_mf_RMSE_xk;...
                   RBPS_dep_rts_RMSE_xk;MHRBPS_dep_rts_RMSE_xk;MHRBPP_dep_rts_RMSE_xk;...
                   KS_dep_RMSE_xk];
Averaged_time=[Timebrbpfdep,Timepsdep,Timerbpsind,...
               Timerbpsdepmf,Timemhrbpsdepmf,Timemhrbppdepmf,...
               Timerbpsdeprts,Timemhrbpsdeprts,Timemhrbppdeprts,...
               Timeksdep]./K;
t=0:m.ss.T-1;
figure(1);
plot(t,BRBPF_dep_RMSE_xp,'-sr',t,PS_dep_RMSE_xp,'-og',t,RBPS_ind_RMSE_xp,'-hm',t,RBPS_dep_mf_RMSE_xp,'-<b',t,MHRBPS_dep_mf_RMSE_xp,'-cd',t,MHRBPP_dep_mf_RMSE_xp,'-+k',t,KS_dep_RMSE_xp,'-py','linewidth',1);
legend('D-RBPF','D-PS','I-RBPS','D-RBPS','D-MH-RBPS1','D-MH-RBPS2','D-KS');
xlabel('Time step');
ylabel('Position RMSE');
figure(2);
plot(t,BRBPF_dep_RMSE_xk,'-sr',t,PS_dep_RMSE_xk,'-og',t,RBPS_ind_RMSE_xk,'-hm',t,RBPS_dep_mf_RMSE_xk,'-<b',t,MHRBPS_dep_mf_RMSE_xk,'-cd',t,MHRBPP_dep_mf_RMSE_xk,'-+k',t,KS_dep_RMSE_xk,'-py','linewidth',1);
legend('D-RBPF','D-PS','I-RBPS','D-RBPS','D-MH-RBPS1','D-MH-RBPS2','D-KS');
xlabel('Time step');
ylabel('Velocity RMSE');

figure(3);
plot(t,RBPS_dep_mf_RMSE_xp,'-<b',t,MHRBPS_dep_mf_RMSE_xp,'-cd',t,MHRBPP_dep_mf_RMSE_xp,'-+k',t,RBPS_dep_rts_RMSE_xp,'-sr',t,MHRBPS_dep_rts_RMSE_xp,'-og',t,MHRBPP_dep_rts_RMSE_xp,'-hm','linewidth',1);
legend('D-RBPS','D-MH-RBPS1','D-MH-RBPS2','D-RBPS*','D-MH-RBPS1*','D-MH-RBPS2*');
xlabel('Time step');
ylabel('Position RMSE');
figure(4);
plot(t,RBPS_dep_mf_RMSE_xk,'-<b',t,MHRBPS_dep_mf_RMSE_xk,'-cd',t,MHRBPP_dep_mf_RMSE_xk,'-+k',t,RBPS_dep_rts_RMSE_xk,'-sr',t,MHRBPS_dep_rts_RMSE_xk,'-og',t,MHRBPP_dep_rts_RMSE_xk,'-hm','linewidth',1);
legend('D-RBPS','D-MH-RBPS1','D-MH-RBPS2','D-RBPS*','D-MH-RBPS1*','D-MH-RBPS2*');
xlabel('Time step');
ylabel('Velocity RMSE');
disp('Pos_ARMSEs for D-RBPF,D-PS,I-RBPS,D-RBPS,D-MH-RBPS1,D-MH-RBPS2,D-RBPS*,D-MH-RBPS1*,D-MH-RBPS2*,D-KS are:')
mean(Position_RMSE_Arr,2)'
disp('Vel_ARMSEs for D-RBPF,D-PS,I-RBPS,D-RBPS,D-MH-RBPS1,D-MH-RBPS2,D-RBPS*,D-MH-RBPS1*,D-MH-RBPS2*,D-KS are:')
mean(Velocity_RMSE_Arr,2)'
disp('A_Run_time for D-RBPF,D-PS,I-RBPS,D-RBPS,D-MH-RBPS1,D-MH-RBPS2,D-RBPS*,D-MH-RBPS1*,D-MH-RBPS2*,D-KS are:')
Averaged_time





