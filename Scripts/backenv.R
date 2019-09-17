h = function(x) {3/4*(1-exp(-4/3*x))}
ed = function(x) {-3/4*log(1-4/3*(x))}

sample <- function (d = 0.2, L = 100, n=1000) {
  hh = rpois(n = n,lambda = h(d)*L)
  #ed(min((hh+0.01)/L,3/4-0.001))
  ed((hh+0.01)/L)
}

require(ggplot2); require(reshape2)

g=expand.grid(d=c(10:20)/100,L=c(100))
g = cbind (g,replicate(4000,apply(g,1,function(x) sample(x[[1]],x[[2]],1))))
g = melt(g,id.vars = 1:2)

qplot(d,value/d,data=g,group=d,geom="violin")+scale_y_log10()+facet_wrap(.~L,ncol = 1)


ggplot(aes(x=x/0.1),data=data.frame(x=sample(d=0.1,L=500,n=10000)))+stat_ecdf()+stat_function(
  fun =function(x,y=0.1,L=500) {ppois(q = round(L*h(x*y)), lambda = h(y)*L )}, color="red")

h=function(x) x
m=100
d=0.1
ll=1000
d*ll
ggplot(data.frame(x = c(1/m, m)), aes(x))+
# stat_function(fun =function(x,y=d,L=ll) {
#    ((-log( dpois(x = round(L*h(x*y) ), lambda = h(y)*L )) +  log( dpois(x = round(L*h(y) ) , lambda = h(y)*L ))))}, color="red")+
#  stat_function(fun =function(x,y=d,L=ll) {
#    ((-log( dpois(x = round(L*h(y) ), lambda = h(y*x)*L )) +  log( dpois(x = round(L*h(y) ) , lambda = h(y)*L ))))}, color="green")+
  stat_function(fun =function(x,y=0.1,L=500) {((log(x)-x)-(log(1)-1))/((log(m)-m)-(log(1)-1))}, color="blue")+
  stat_function(fun =function(x,y=0.1,L=500) {log(x)*log(x)/(log(m)*log(m))}, color="black")
#stat_function(fun =function(x,y=0.1,L=500) {(1/d*(1-x)^2)/(1/d*(1-m)^2)}, color="blue")

c(1/2*qchisq(p=0.05,df=2*ll*d)/(ll*d),1/2*qchisq(p=0.95,df=2*ll*d+2)/(ll*d))

ggplot(data.frame(x = c(0.001,2)), aes(x))+geom_hline(yintercept = 1,color="red")+
  stat_function(fun=function(x) 1/2*qchisq(p=0.05,df=2*ll*x)/(ll*x))+
  stat_function(fun=function(x) 1/2*qchisq(p=0.95,df=2*ll*x+2)/(ll*x))+
  stat_function(fun=function(x) 1/2*qchisq(p=0.2,df=2*ll*x)/(ll*x))+
  stat_function(fun=function(x) 1/2*qchisq(p=0.8,df=2*ll*x+2)/(ll*x))+
  scale_y_log10()+scale_x_log10()
