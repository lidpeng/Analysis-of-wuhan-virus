infect = t(read_csv("confirmed_cases.csv"))
cor.test()
plot(Total,type="l",main="总感染人数",xlab="Time")
dif = diff(Total)
plot(dif,type="l",main="每个时间点新增感染人数",xlab="Time")
dif2 = diff(dif)
plot(dif2,type="l",main="总感染人数的二阶差分图",xlab="Time")
for(i in 1:4)print(Box.test(dif2,type = "Ljung-Box",lag = 6*i)) ## 白噪声检验

a1=arima(dif2,order=c(0,1,1))
for(i in 1:4)print(Box.test(a1$residuals,lag = 6*i,type = "Ljung-Box"))
a2=arima(dif2,order=c(5,1,0))
for(i in 1:4)print(Box.test(a2$residuals,lag = 6*i,type = "Ljung-Box"))
a3=arima(dif2,order=c(5,1,1))
for(i in 1:4)print(Box.test(a3$residuals,lag = 6*i,type = "Ljung-Box"))

real=c(14525,14536,16639,16640,16651,16667,16676,16691,16832,16901,16917,17147,17170,17184,17187,17200)
difreal = diff(diff(real))##截止2-3 10:56
predict(a2,14)

getFitData <- function(x,
                       cond.dist = "QMLE",
                       include.mean = FALSE,
                       ...)
{
  args <- list(...)
  args$data = x
  args$cond.dist = cond.dist
  args$include.mean = include.mean
  
  log <- capture.output(
    {
      fit <- do.call(garchFit, args = args)
    })
  
  res <- coef(fit)
  res[paste0(names(fit@fit$se.coef), ".se")] <- fit@fit$se.coef
  return(res)
}
experiments <- foreach(
  s = seeds) %do%
  {
    set.seed(s)
    x <- garchSim(
      garchSpec(
        model = list(
          alpha = 0.2, beta = 0.2, omega = 0.2)),
      n.start = 1000, n = 1000)
    params <- foreach(
      t = 50:1000,
      .combine = rbind,
      .packages = c("fGarch")) %dopar%
      {
        getFitData(x[1:t])
      }
    rownames(params) <- 50:1000
    params
  }
names(experiments) <- seeds

q=seq(from=1,to=11,by=1)
plot(daily~q)
lm1=lm(daily~q)
lm2=lm(daily~q+I(q^2))
lm3=lm(daily~q+I(q^2)+I(q^3))
lm4=lm(daily~q+I(q^2)+I(q^3)+I(q^4))
lines(q,fitted(lm1),col='green',lwd=2)
lines(q,fitted(lm2),col='red',lwd=2)
lines(q,fitted(lm3),col='black',lwd=2)
lines(q,fitted(lm4),col='orange',lwd=2)
legend("bottomright",c("Linear","Quadratic","Cubic","Quartic"),
       col=c("green","red","black","orange"), lwd=3)

def1=abs(fitted(lm1)-daily)/daily
def2=abs(fitted(lm2)-daily)/daily
def3=abs(fitted(lm3)-daily)/daily
def4=abs(fitted(lm4)-daily)/daily
plot(q,def1,type="l",col='green')
lines(q,def2,col='red',lwd=2)
lines(q,def3,col='black',lwd=2)
lines(q,def4,col='orange',lwd=2)

lmf1=function(x){
  y=lm1$coefficients[1]+lm1$coefficients[2]*x
  return(y)
}
lmf2=function(x){
  y=lm2$coefficients[1]+lm2$coefficients[2]*x+lm2$coefficients[3]*x^2
  return(y)
}
lmf3=function(x){
  y=lm3$coefficients[1]+lm3$coefficients[2]*x+lm3$coefficients[3]*x^2+lm3$coefficients[4]*x^3
  return(y)
}
lmf4=function(x){
  y=lm4$coefficients[1]+lm4$coefficients[2]*x+lm4$coefficients[3]*x^2+lm4$coefficients[4]*x^3+lm4$coefficients[5]*x^4
  return(y)
}
plot(seq(1,18,1),daily,ylim=c(0,36200),type='b')
q=seq(1,18,1)
lines(seq(11,18,1),lmf1(11:18),col='green',lwd=2)
lines(seq(11,18,1),lmf2(11:18),col='red',lwd=2)
lines(seq(11,18,1),lmf3(11:18),col='black',lwd=2)
lines(seq(11,18,1),lmf4(11:18),col='orange',lwd=2)
legend("topleft",c("Linear","Quadratic","Cubic","Quartic"),
       col=c("green","red","black","orange"), lwd=3)

coo = read.csv("cor.csv",header=T,row.names=1)
ggcorrplot(coo, hc.order = TRUE, outline.color = "white") 