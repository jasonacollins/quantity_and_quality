#This code runs the three genotype simulation contained in Collins, Jason, Boris Baer and E. Juerg Weber (2013) Economic Growth and Evolution: Parental Preferences for Quality and Quantity of Offspring, Macroeconomic Dynamics (forthcoming)

# For further information, contact Jason Collins at jason@jasoncollins.org, or visit Evolving Economics http://www.jasoncollins.org

#load BB package - used for solving nonlinear equations
library(BB)

#initial conditions
A<-1
g<-0
ea<-0
eb<-0
ec<-0
e<-0
La<-0.007
Lb<-0.7
Lc<-0.007
L<-La+Lb+Lc
qa<-La/L
qb<-Lb/L
qc<-Lc/L
za<-1.25
zb<-1.25
zc<-1.25
na<-1
nb<-1
nc<-1

#parameters
Ba<-1
Bb<-0.9
Bc<-0.75
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
Growth<-data.frame(time=0, A, g, ea, eb, ec, La, Lb, Lc, na, nb, nc, za, zb, zc)

#establish a loop
for (t in 1:time) {

# population
La<-na*La
Lb<-nb*Lb
Lc<-nc*Lc
L<-La+Lb+Lc
qa<-La/L
qb<-Lb/L
qc<-Lc/L

# technology in this period (based on education given to children in last period)
e<-qa*ea+qb*eb+qc*ec
g<-k*e^0.5
A<-(1+g)*A

# human capital
ha<-(m*ea+a)/(ea+r*g+a)
hb<-(m*eb+a)/(eb+r*g+a)
hc<-(m*ec+a)/(ec+r*g+a)

#level of eduction of each genotype
ea<-max(0, (1/(2*m))*((Ba*m*r*g+Ba*m*a-Ba*a-m*r*g-a*m-a)+(((Ba*m*r*g+Ba*m*a-Ba*a-m*r*g-a*m-a)^2+4*m*(Ba*m*r*g*Tau+Ba*m*a*Tau-Ba*a*Tau-a*r*g-a^2))^0.5)))
eb<-max(0, (1/(2*m))*((Bb*m*r*g+Bb*m*a-Bb*a-m*r*g-a*m-a)+(((Bb*m*r*g+Bb*m*a-Bb*a-m*r*g-a*m-a)^2+4*m*(Bb*m*r*g*Tau+Bb*m*a*Tau-Bb*a*Tau-a*r*g-a^2))^0.5)))
ec<-max(0, (1/(2*m))*((Bc*m*r*g+Bc*m*a-Bc*a-m*r*g-a*m-a)+(((Bc*m*r*g+Bc*m*a-Bc*a-m*r*g-a*m-a)^2+4*m*(Bc*m*r*g*Tau+Bc*m*a*Tau-Bc*a*Tau-a*r*g-a^2))^0.5)))

#calculate possibilities for income for genotype a
#if subsistence constraint binding for noone
zaa1b1c1<-(((A*X)/(L*(1-gamma)*(qa*ha+qb*hb+qc*hc)))^alpha)*ha
zba1b1c1<-zaa1b1c1*hb/ha
zca1b1c1<-zaa1b1c1*hc/ha
if(zaa1b1c1>=zsc & zba1b1c1>=zsc & zca1b1c1>=zsc) za<-zaa1b1c1

#if subsistence constraint binding for type c only
fzaa1b1c2<-function(zaa1b1c2){
y<-(((A*X)/(L*(((1-gamma)*(qa*ha+qb*hb))+qc*((sc*ha)/zaa1b1c2))))^alpha)*ha-zaa1b1c2
y
}
ansfzaa1b1c2<-multiStart(c(0.2,1,10,100),fzaa1b1c2, control=list(M=200), quiet=TRUE)
zaa1b1c2<-max(ansfzaa1b1c2$par[,1], ansfzaa1b1c2$par[,2], ansfzaa1b1c2$par[,3], ansfzaa1b1c2$par[,4])
zba1b1c2<-zaa1b1c2*hb/ha
zca1b1c2<-zaa1b1c2*hc/ha
if(zaa1b1c2>=zsc & zba1b1c2>=zsc & zca1b1c2>=sc & zca1b1c2<zsc) za<-zaa1b1c2

#if subsistence constraint binding for type b and c
fzaa1b2c2<-function(zaa1b2c2){
y<-(((A*X)/(L*(((1-gamma)*qa*ha)+qb*((sc*ha)/zaa1b2c2)+qc*((sc*ha)/zaa1b2c2))))^alpha)*ha-zaa1b2c2
y
}
ansfzaa1b2c2<-multiStart(c(0.2,1,10,100),fzaa1b2c2, control=list(M=200), quiet=TRUE)
zaa1b2c2<-max(ansfzaa1b2c2$par[,1], ansfzaa1b2c2$par[,2], ansfzaa1b2c2$par[,3], ansfzaa1b2c2$par[,4])
zba1b2c2<-zaa1b2c2*hb/ha
zca1b2c2<-zaa1b2c2*hc/ha
if(zaa1b2c2>=zsc & zba1b2c2>=sc & zba1b2c2<zsc & zca1b2c2>=sc & zca1b2c2<zsc) za<-zaa1b2c2

#if subsistence constraint binding for all types
fzaa2b2c2<-function(zaa2b2c2){
y<-(((A*X)/(L*(((sc/zaa2b2c2)*qa*ha)+(1-qa)*((sc*ha)/zaa2b2c2))))^alpha)*ha-zaa2b2c2
y
}
ansfzaa2b2c2<-multiStart(c(0.2,1,10,100),fzaa2b2c2, control=list(M=200), quiet=TRUE)
zaa2b2c2<-max(ansfzaa2b2c2$par[,1], ansfzaa2b2c2$par[,2], ansfzaa2b2c2$par[,3], ansfzaa2b2c2$par[,4])
zba2b2c2<-zaa2b2c2*hb/ha
zca2b2c2<-zaa2b2c2*hc/ha
if(zaa2b2c2>=sc & zaa2b2c2<zsc & zba2b2c2>=sc & zba2b2c2<zsc) za<-zaa2b2c2
if(zaa2b2c2<sc) za<-0 

#given income for genotype a, calculate income for genotype b
zb<-za*hb/ha
zc<-za*hc/ha

#population growth
if(za>=zsc) na<-gamma/(Tau+ea)
if(za<zsc & za>sc) na<-(1-(sc/za))/(Tau+ea)
if(za<=sc) na<-0
if(zb>=zsc) nb<-gamma/(Tau+eb)
if(zb<zsc & zb>sc) nb<-(1-(sc/zb))/(Tau+eb)
if(zb<=sc) nb<-0
if(zc>=zsc) nc<-gamma/(Tau+ec)
if(zc<zsc & zc>sc) nc<-(1-(sc/zc))/(Tau+ec)
if(zc<=sc) nc<-0

Growth<-rbind(Growth, c(t, A, g, ea, eb, ec, La, Lb, Lc, na, nb, nc, za, zb, zc))

}

write.table(Growth,file="Growth.txt", sep = ",")

Growth