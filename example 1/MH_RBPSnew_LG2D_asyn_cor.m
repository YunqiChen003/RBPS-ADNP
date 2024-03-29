function G=MH_RBPSnew_LG2D_asyn_cor(M,m)
% M: output of Rao-Blackwellzed particle filters;
% m: model parametes;
fprintf('RBPS_asyn_cor\n');
G1=MH_RBBP_LG2Dasyn_cor1(M,m);
%% MH_RBPSnew using conditional MF smoothing;
G2=Conditional_MF_smoothing_asyn_cor(G1,m);
G.time1=G1.time+G2.time;
G.xhat1=[G1.xnhat;G2.xlhat];
G.Xerror1=m.x-G.xhat1;
%% MH_RBPSnew using conditional RTS smoothing;
G3=Conditional_RTS_smoothing_asyn_cor1(G1,m);
G.time2=G1.time+G3.time;
G.xhat2=[G1.xnhat;G3.xlhat];
G.Xerror2=m.x-G.xhat2;


