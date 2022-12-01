endT=50;tstep=1
pop0=c(990,10,0,0) #SIRV
pars=list(c=4,p=0.1,ve=0.2,d=5)
n=100 #number of iterations
costs=list(I=100,vx=10,q=400,s=100)
library(tidyverse)


## Tau Leap
Rates_tau=function(params,x,tstep){
  S=x[1];I=x[2];R=x[3];V=x[4];N=S+I+R+V
  rates=rep(NA,2)
  rates[1]=params$beta*S*I/N #infection
  rates[2]=params$gamma*I #recovery
  rates[3]=params$beta_v*S*V/N #vaccinated infection
  rates=rates*tstep #multiply all rates by time step
  return(rates)
}
Events_tau=function(rates,x){
  S=x[1];I=x[2];R=x[3];V=x[4];N=S+I+R+V
  infection=min(S,rpois(1,rates[1])) #we don't want more infected than the total in S
  S=S-infection;I=I+infection #infection
  infection_v=min(V,rpois(1,rates[3])) #we don't want more vx infected than the total in V
  V=V-infection_v;I=I+infection_v #vx infection
  recovery=min(I,rpois(1,rates[2])) #we don't want more recovering than the total in I
  I=I-recovery;R=R+recovery #recovery
  return(c(S,I,R,V))
}
Iteration_tau=function(pars,tstep,endT,pop0,q,v,s){
  pop_tau=matrix(NA,nrow=endT/tstep+1,ncol=5,dimnames = list(NULL,c("time","S","I","R","V")))
  vx=rep(0,nrow(pop_tau))
  pop_now=pop0
  t=0;step=1
  pop_tau[1,]=c(t,pop_now)
  while(t<endT){
    if(v[t+1]>0){if(v[t+1]<pop_now[1]){
      vx[t]=v[t];pop_now[1]=pop_now[1]-v[t+1];pop_now[4]=pop_now[4]+v[t+1]}else{
        vx[t]=pop_now[1];pop_now[1]=0;pop_now[4]=pop_now[4]+vx[t+1]
      }}
    params=list(beta=pars$c*pars$p*q[t+1],beta_v=pars$c*pars$p*pars$ve*q[t+1],gamma=s[t+1]/pars$d) #parameters
    ratedraw=Rates_tau(params,pop_now,tstep) #calculate rates
    eventdraw=Events_tau(ratedraw,pop_now) #make events happen
    t=t+tstep;step=step+1 #update time and matrix position
    pop_tau[step,]=c(t,eventdraw)
    pop_now=eventdraw #update population
  }
  return(list(finalpop=pop_now,v_n=sum(vx),I_n=sum(pop_tau[,3])-1,q_n=length(which(q!=1)),s_n=length(which(s!=1))))
}

Cost_est=function(endT,pop0,pars,n,q,v,s,costs){
pop_sim_tau=list(NULL);
outcomes=data.frame(total_infection=rep(NA,n),total_vx=rep(NA,n),total_q=rep(NA,n),total_s=rep(NA,n))
for(i in 1:n){
  pop_sim_tau[[i]]=Iteration_tau(pars,tstep,endT,pop0,q,v,s) #results of iteration
  outcomes$total_infection[i]=pop_sim_tau[[i]]$I_n
  outcomes$total_vx[i]=pop_sim_tau[[i]]$v_n
  outcomes$total_q[i]=pop_sim_tau[[i]]$q_n
  outcomes$total_s[i]=pop_sim_tau[[i]]$s_n
}
outcomes$costs=outcomes$total_infection*costs$I+outcomes$total_vx*costs$vx+outcomes$total_q*costs$q+outcomes$total_s*costs$s
return(outcomes)
}


#Basic scenarios
qs=cbind(rep(1,endT),rep(0.5,endT),rep(1,endT),rep(1,endT),rep(0.5,endT))
vs=cbind(rep(0,endT),rep(0,endT),rep(20,endT),rep(0,endT),rep(20,endT))
ss=cbind(rep(1,endT),rep(1,endT),rep(1,endT),rep(2,endT),rep(2,endT))
out_cost=matrix(NA,nrow=n,ncol=ncol(qs))
colnames(out_cost)=c("Null","Q","V","S","QVS")

#Class scenarios
qs=cbind(c(rep(0.5,7),rep(1,endT-7)),
         c(rep(1,20),rep(0.5,5),rep(1,15),rep(0.5,10)),
         c(rep(0.5,30),rep(1,20)),
         rep(1,endT),
         c(rep(0.5,5),rep(1,45)))
vs=cbind(rep(20,endT),
         c(rep(20,20),rep(0,endT-20)),
         rep(20,endT),
         rep(20,endT),
         c(rep(0,5),rep(20,25),rep(0,20)))
ss=cbind(c(rep(1,5),rep(2,20),rep(1,endT-25)),
         c(rep(1,10),rep(2,15),rep(1,endT-25)),
         c(rep(1,40),rep(2,10)),
         rep(1,endT),
         c(rep(1,30),rep(2,20)))
out_cost=matrix(NA,nrow=n,ncol=ncol(qs))
colnames(out_cost)=c("E",
                     "A",
                     "S",
                     "M",
                     "T")

costs=list(I=50,vx=10,q=200,s=100)

for(p in 1:ncol(qs)){
  out=Cost_est(endT,pop0,pars,n,qs[,p],vs[,p],ss[,p],costs)
  out_cost[,p]=out$costs
}
cost_long=pivot_longer(as.data.frame(out_cost),names_to = "Scenario",values_to = "Cost",cols = 1:ncol(out_cost))
ggplot(cost_long,aes(x=Scenario,y=Cost))+geom_boxplot()+geom_jitter()

