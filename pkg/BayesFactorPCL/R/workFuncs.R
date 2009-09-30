

################
# Work Functions
################

jnt.post.norm = function(eta2, delta, t, N, r.scale=1)
exp(
-(N/2+1)*log(eta2) - (N-1)/(2*eta2) - (t - delta*sqrt(eta2*N))^2/(2*eta2) - delta^2/(2*r.scale^2)
)


jnt.post = function(eta2, delta, t, N, df.prior=1, r.scale=1)
exp(
-(N/2+1)*log(eta2) - (N-1)/(2*eta2) - (t - delta*sqrt(eta2*N))^2/(2*eta2) - (df.prior+1)/2*log(1+(delta/r.scale)^2/df.prior)
)

vect.jnt.post = function(eta2, delta, t, N, df.prior=1,r.scale=1){
if(df.prior==Inf){
return(sapply(eta2,jnt.post.norm,delta=delta,t=t,N=N,r.scale=r.scale))
}else{
return(sapply(eta2,jnt.post,delta=delta,t=t,N=N,df.prior=df.prior,r.scale=r.scale))
}
}




vect.area.bf = Vectorize(area.bf,c("t","N1"))


margd.post = function(delta, t, N, df.prior=1,eta2.lo=0,eta2.up=Inf,r.scale=1)
integrate(vect.jnt.post, eta2.lo, eta2.up, delta=delta,t=t,N=N,df.prior=df.prior,r.scale=r.scale,rel.tol=max(50*.Machine$double.eps,0.5e-28))$value

vect.margd.post = function(delta, t, N, df.prior=1,eta2.lo=0,eta2.up=Inf,r.scale=1)
sapply(delta, margd.post, t=t,N=N,df.prior=df.prior,eta2.lo=eta2.lo,eta2.up=eta2.up,r.scale=r.scale)


critt.work=function(t,N,null.d=0.2,df.prior=1,bf=1/10,r.scale=1)
(log(bf)-log(vect.calc.bf.area2(t,N,null.d=0.2,df.prior=1,r.scale=r.scale)))^2


find.crit.t.areabf = function(bf=1/10,N,null.d=0.2,df.prior=1,bounds=c(0,20),r.scale=1)
optimize(critt.work,interval=bounds,N=N,null.d=null.d,df.prior=df.prior,r.scale=r.scale)$minimum

vect.crit.t.areabf=Vectorize(find.crit.t.areabf,c("N"))


sig2ciUp=function(s2,N,alpha=0.05)
s2*(N-1)/qchisq(alpha,N-1)


s2.interval=function(y,alpha=.05){
N=length(y)
s2=var(y)
cutoffhi=qchisq(alpha/2,N-1)
CIhi=s2*(N-1)/cutoffhi
cutofflo=qchisq(1-alpha/2,N-1)
CIlo=s2*(N-1)/cutofflo
c(CIlo,CIhi)
}

proplog.ddelta.t=function(d,eta2,t,N,nu=1,r=1)
-(nu+1)/2*log(1+(d/r)^2/nu)-1/(2*eta2)*(t-d*sqrt(N*eta2))^2


proplog.deta2.t=function(eta2,d,t,N)
(-N/2-1)*log(eta2)-1/(2*eta2)*(t-d*sqrt(N*eta2))^2 - (N-1)/(2*eta2)


samp.delta.t=function(old,eta2,t,N,sdmet,nu=1,r=1){
        if(nu==Inf){
                return(rnorm(1,sqrt(N)/(N+1/(r^2))*t/sqrt(eta2),sqrt(1/(N+1/(r^2)))))
        }else{
                z=rnorm(1,0,sdmet)
                cand=old+z
                llr=proplog.ddelta.t(cand,eta2,t,N)-proplog.ddelta.t(old,eta2,t,N,r=r)
                b=rbinom(1,1,exp(min(llr,0)))
                return(b*cand+(1-b)*old)
        }
}

samp.eta2.t=function(old,d,t,N,sdmet){
z=rnorm(1,0,sdmet)
cand=old+z
if(cand<=0){
llr=-Inf
}else{
llr=proplog.deta2.t(cand,d,t,N)-proplog.deta2.t(old,d,t,N)
}
b=rbinom(1,1,exp(min(llr,0)))
b*cand+(1-b)*old
}


find.ncp.lo=function(ncp,alpha,t,N)
(t-qt(1-alpha/2,N-1,ncp=ncp))^2

find.ncp.hi=function(ncp,alpha,t,N)
(t-qt(alpha/2,N-1,ncp=ncp))^2

freqCIdelta=function(t,N,alpha=.05){
ow=options("warn")
options(warn=-1)
t.hi=optimize(find.ncp.hi,interval=c(-10,10),tol = .Machine$double.eps^0.25,t=t,alpha=alpha,N=N)$minimum
t.lo=optimize(find.ncp.lo,interval=c(-10,10),tol = .Machine$double.eps^0.25,t=t,alpha=alpha,N=N)$minimum

d.hi=t.hi/sqrt(N)
d.est=t/sqrt(N)
d.lo=t.lo/sqrt(N)

options(ow)
c(d.lo,d.est,d.hi)
}


