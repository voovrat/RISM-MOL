function Scheme = GenerateNMGMScheme(eps, mu, SmoothSteps, CoarseSteps )

INF_ITER = 100000;

Scheme = [ struct('N',4096,'dR',0.05,'Niter_pre',SmoothSteps,'Niter_post',SmoothSteps,'eps',eps,'mu',1)
           struct('N',2048,'dR',0.05,'Niter_pre',SmoothSteps,'Niter_post',SmoothSteps,'eps',eps,'mu',INF_ITER)
           struct('N',1024,'dR',0.05,'Niter_pre',SmoothSteps,'Niter_post',SmoothSteps,'eps',eps,'mu',mu)
           struct('N',512,'dR',0.1,'Niter_pre',SmoothSteps,'Niter_post',SmoothSteps,'eps',eps,'mu',mu)
           struct('N',256,'dR',0.2,'Niter_pre',SmoothSteps,'Niter_post',SmoothSteps,'eps',eps,'mu',mu) 
           struct('N',128,'dR',0.4,'Niter_pre',CoarseSteps,'Niter_post',0,'eps',eps,'mu',0) ];
