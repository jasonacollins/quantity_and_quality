#This code runs the two genotype simulation contained in Collins, Jason, Boris Baer and E. Juerg Weber (2013) Economic Growth and Evolution: Parental Preferences for Quality and Quantity of Offspring, Macroeconomic Dynamics (forthcoming)

# For further information, contact Jason Collins at jason@jasoncollins.org, or visit Evolving Economics http://www.jasoncollins.org

#load BB package - used for solving nonlinear equations
library(BB)

#initial conditions
A<-1
g<-0
ea<-0
eb<-0
e<-0
La<-0.007
Lb<-0.7
L<-La+Lb
qa<-La/L
qb<-Lb/L
za<-1.25
zb<-1.25
na<-1
nb<-1

#parameters
Ba<-1
Bb<-0.9
alpha<-0.4
Tau<-0.2
rho<-0.99
a<-rho*Tau
m<-2
gamma<-0.259
k<-8.885139596
r<-0.108150721
X<-1
sc<-1
zsc<-sc/(1-gamma)
time<-200 #number of generations

#Build data frame which will be used to store results
Growth<-data.frame(time=0, A, g, ea, eb, La, Lb, na, nb, za, zb)

#establish a loop
for (t in 1:time) {

# population
La<-na*La
Lb<-nb*Lb
L<-La+Lb
qa<-La/L
qb<-Lb/L

# technology in this period (based on education given to children in last period)
e<-qa*ea+qb*eb
g<-k*e^0.5
A<-(1+g)*A

# human capital
ha<-(m*ea+a)/(ea+r*g+a)
hb<-(m*eb+a)/(eb+r*g+a)

#level of eduction of each genotype
ea<-max(0, (1/(2*m))*((Ba*m*r*g+Ba*m*a-Ba*a-m*r*g-a*m-a)+(((Ba*m*r*g+Ba*m*a-Ba*a-m*r*g-a*m-a)^2+4*m*(Ba*m*r*g*Tau+Ba*m*a*Tau-Ba*a*Tau-a*r*g-a^2))^0.5)))
eb<-max(0, (1/(2*m))*((Bb*m*r*g+Bb*m*a-Bb*a-m*r*g-a*m-a)+(((Bb*m*r*g+Bb*m*a-Bb*a-m*r*g-a*m-a)^2+4*m*(Bb*m*r*g*Tau+Bb*m*a*Tau-Bb*a*Tau-a*r*g-a^2))^0.5)))

#calculate possibilities for income for genotype a (and given zb=za*hb/ha)
#if subsistence constraint binding for neither
zaa1b1<-(((A*X)/(L*(1-gamma)*(qa*ha+(1-qa)*hb)))^alpha)*ha
zba1b1<-zaa1b1*hb/ha
if(zaa1b1>=zsc & zba1b1>=zsc) za<-zaa1b1

#if subsistence constraint binding for type b
fzaa1b2<-function(zaa1b2){
y<-(((A*X)/(L*(((1-gamma)*qa*ha)+(1-qa)*sc/zaa1b2*ha)))^alpha)*ha-zaa1b2
y
}
ansfzaa1b2<-multiStart(c(0.2,1,10,100),fzaa1b2, control=list(M=200), quiet=TRUE)
zaa1b2<-max(ansfzaa1b2$par[,1], ansfzaa1b2$par[,2], ansfzaa1b2$par[,3], ansfzaa1b2$par[,4])
zba1b2<-zaa1b2*hb/ha
if(zaa1b2>=zsc & zba1b2>=sc & zba1b2<zsc) za<-zaa1b2

#if subsistence constraint binding for both types
fzaa2b2<-function(zaa2b2){
y<-(((A*X)/(L*((sc/zaa2b2*qa*ha)+(1-qa)*sc/zaa2b2*ha)))^alpha)*ha-zaa2b2
y
}
ansfzaa2b2<-multiStart(c(0.2,1,10,100),fzaa2b2, control=list(M=200), quiet=TRUE)
zaa2b2<-max(ansfzaa2b2$par[,1], ansfzaa2b2$par[,2], ansfzaa2b2$par[,3], ansfzaa2b2$par[,4])
zba2b2<-zaa2b2*hb/ha
if(zaa2b2>=sc & zaa2b2<zsc & zba2b2>=sc & zba2b2<zsc) za<-zaa2b2
if(zaa2b2<sc) za<-0 

#given income for genotype a, calculate income for genotype b
zb<-za*hb/ha

#population growth
if(za>=zsc) na<-gamma/(Tau+ea)
if(za<zsc & za>sc) na<-(1-(sc/za))/(Tau+ea)
if(za<=sc) na<-0
if(zb>=zsc) nb<-gamma/(Tau+eb)
if(zb<zsc & za>sc) nb<-(1-(sc/zb))/(Tau+eb)
if(zb<=sc) nb<-0

Growth<-rbind(Growth, c(t, A, g, ea, eb, La, Lb, na, nb, za, zb))

}

write.table(Growth,file="Growth.txt", sep = ",")

Growth