function e=MHproposal(xparticle1,xlfilter1,Plfilter1,xparticle2,m)
mu1=xparticle1+m.ss.A(1,2)*xlfilter1;
P11=m.ss.An*Plfilter1*m.ss.An'+m.ss.Qnn;
mu2=m.ss.Al*xlfilter1;
P22=m.ss.Al*Plfilter1*m.ss.Al'+m.ss.Qll;
P12=m.ss.An*Plfilter1*m.ss.Al'+m.ss.Qnl;
Cov_13=P11+P12*m.ss.An';
P33=(P11+m.ss.An*P12')+(P12+m.ss.An*P22)*m.ss.An'+m.ss.Qnn;
mu3=mu1+m.ss.An*mu2;
K=Cov_13/P33;
mu1_3=mu1+K*(xparticle2-mu3);
P11_3=P11-K*P33*K';
e.mu=mu1_3;
e.P=P11_3;



