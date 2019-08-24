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

m=2
d=0.1
ll=500
ggplot(data.frame(x = c(1/m, m)), aes(x))+stat_function(fun =function(x,y=d,L=ll) {(-log( dpois(x = round(L*h(x*y) ), lambda = h(y)*L )) +  log( dpois(x = round(L*h(y) ) , lambda = h(y)*L )))/(-log( dpois(x = round(L*h(m*y) ), lambda = h(y)*L )) +  log( dpois(x = round(L*h(y) ) , lambda = h(y)*L ))) }, color="red")+
  stat_function(fun =function(x,y=0.1,L=500) {(1/d*(1-x)^2)/(1/d*(1-m)^2)}, color="blue")+
  stat_function(fun =function(x,y=0.1,L=500) {log(x)*log(x)/(log(m)*log(m))}, color="black")



