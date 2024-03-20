function G=RBPS_asyn_cor(M,m)
% M: output of Rao-Blackwellzed particle filters;
% m: model parametes;
fprintf('RBPS_asyn_cor\n');
G1=RB_backward_simluation_asyn_cor(M,m);
%% RBPS using conditional MF smoothing;
G2=Conditional_MF_smoothing_asyn_cor2(G1,m);
G.time1=G1.time+G2.time;
G.xhat1=[G1.xnhat;G2.xlhat];
G.Xerror1=m.x-G.xhat1;
%% RBPS using conditional RTS smoothing;
G3=Conditional_RTS_smoothing_asyn_cor2(G1,m);
G.time2=G1.time+G3.time;
G.xhat2=[G1.xnhat;G3.xlhat];
G.Xerror2=m.x-G.xhat2;



