function G=RBPS_ind(M,m)
% M: output of Rao-Blackwellzed particle filters;
% m: model parametes;
fprintf('RBPS_ind\n');
G1=RB_backward_simluation_ind(M,m);
G2=Conditional_MF_smoothing_ind2(G1,m);
G.time=G1.time+G2.time;
G.xhat=[G1.xnhat;G2.xlhat];
G.Xerror=m.x-G.xhat;
