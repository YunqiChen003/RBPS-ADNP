function e=trans_likeli(xparticle,xparticle1,xlfilter1,Plfilter1,m,t)
% calculate  p(x_{k}^n | x_{0:k-1}^n, y_{0:k-1}) and p(y_{k} | x_{0:k}^n, y_{0:k-1}), 
% xparticle: dimx*m.ss.Ns particles of x_{k}^n;
% xparticle1：dimx*m.ss.Ns particles of x_{k-1}^{n};
% xfilter：dimx*m.ss.Ns means of x_{k-1}^{l};
% Plfilter1: dimxl*dimxl covariance of x_{k-1}^{l};
% m: model parameters;
% t: time step;
%% PF time update
e.mux=xparticle1+m.ss.An*xlfilter1;
e.Sigmaxx=m.ss.An*Plfilter1*m.ss.An'+m.ss.Qnn;
e.tranP=loggausspdf2(xparticle-e.mux,inv(e.Sigmaxx));
%% KF dynamical measurement update;
    L=Plfilter1*m.ss.An'/e.Sigmaxx;    
    Xlfilterstar=xlfilter1+L*(xparticle-e.mux);
    Plfilterstar=Plfilter1-L*e.Sigmaxx*L';
%% KF time update;
    flbar=m.ss.gamanl'*(xparticle-xparticle1);
    Xlpre=flbar+m.ss.Albar*Xlfilterstar;
    Plpre=m.ss.Albar*Plfilterstar*m.ss.Albar'+m.ss.Qllbar;
%% KF measurement update
    Skesai=[Plfilterstar,Plfilterstar*m.ss.Albar';m.ss.Albar*Plfilterstar,Plpre];
    Sigmayy=m.ss.Hbar*Skesai*m.ss.Hbar'+m.ss.R0;%PF measurement update;
    muy=mea_eq(xparticle)+m.ss.gamany'*(xparticle-xparticle1)-m.ss.gamalybar'*flbar+m.ss.Hbar*[Xlfilterstar;Xlpre];
    K=Skesai(2,:)*m.ss.Hbar'/Sigmayy;%KF measurement update;
    innovation=m.y(:,t)-muy;
    e.xlfilter=Xlpre+K*(innovation);
    e.Plfilter=Plpre-K*Sigmayy*K';
    e.likeP=loggausspdf2(innovation,inv(Sigmayy));