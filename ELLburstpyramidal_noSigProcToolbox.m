%Model of bursting ELL pyramidal cell published: "Dendritic na+ current inactivation can increase cell excitability by delaying a somatic depolarizing afterpotential.
%J Neurophysiol. 2005 Dec;94(6):3836-48." author Fernando R. Fernandez 

clear;
C=1.2; %orig 0.01
Cb=3.5; %orig 0.025
s = 1;
dt=.005;
tt=150
t=0:dt:tt;
I = zeros(s,length(t));
k=1:1:length(t);
kap=.35;
g=1.5;%1.5;

gKmax=10;%15
gNamax =60;%orig 60
gleak=0.18; %orig 0.18

gKmaxb=0; %orig 0
gNamaxb =20;%orig 20
gleakb=0.18;%orig 0.18
gKv3bmax = 8; %orig 8


init =14.5;
step =.5;

% init = .336;

% step = 0.005;

Ek = -88.5;
ENa = 40;
Eleak = -72;

Vc = zeros(s,length(t));
Vc(:,:) = -65;
dVdt= zeros(s,length(t));
Vcb = zeros(s,length(t));
Vcb(:,:) = -65;


taun = zeros(s,length(t));

taum = zeros(s,length(t));

tauh = zeros(s,length(t));

tauk = zeros(s,length(t));


taukb = zeros(s,length(t));

taunb = zeros(s,length(t));

taumb = zeros(s,length(t));

tauhb = zeros(s,length(t));


ninf=zeros(s,length(t));
minf=zeros(s,length(t));
hinf=ones(s,length(t));
kinf=zeros(s,length(t));

nbinf=zeros(s,length(t));
mbinf=zeros(s,length(t));
hbinf=ones(s,length(t));
kbinf=zeros(s,length(t));

n = zeros(s,length(t));
m = zeros(s,length(t));
h = ones(s,length(t));
h(:,:) = .17;
k = zeros(s,length(t));

nb = zeros(s,length(t));
mb = zeros(s,length(t));
hb = ones(s,length(t));
h(:,:) = .17;
kb = zeros(s,length(t));

Ids = zeros(s,length(t));
IKv3 = zeros(s,length(t));
IK = zeros(s,length(t));
INa = zeros(s,length(t));
Ileak = zeros(s,length(t));
Ileak(:,:) = gleak*(-70-Eleak);
Ip = zeros(s,length(t));

Isd = zeros(s,length(t));
IKb = zeros(s,length(t));
INab = zeros(s,length(t));
Ileakb = zeros(s,length(t));
Ileakb(:,:) = gleak*(-70-Eleak);
IKv3b = zeros(s,length(t));

ns=randn(1,length(t));
%[b a]=butter(8,0.02565);%[b a]=butter(8,0.024);
%nsb=filtfilt(b,a,ns);


%----------------------------------------------------------------SOMA

for j = 1:1:s;
   h(:,:) = .12;
   hb(:,:) = .12;
% I(j,1:length(k)) = init + step*j;

count = 1;
spike=1;
threshold =1;
trough=1;
peak =1;
diffwave= 1:1:(4./dt);
for i = 2:1:length(t)-1;
%   I(j,i) = (init + step*j); 
  I(j,i) = (init + step*(j-1)); %+((1*nsb(i)))*0;  

         m(j,i+1) = (1/(1+exp(-(Vc(j,i)+40)/3)));

        
        tauh(j,i+1) = ((2*232/pi)*(28/(4*(Vc(j,i)+64)^2 + 28^2)));

        hinf(j,i+1) = (1/(1+exp(-(Vc(j,i)+40)/-3)));

        h(j,i+1) = (hinf(j,i)-((hinf(j,i)-h(j,i)).*(exp(-(dt)./(tauh(j,i))))));
        
         INa(j,i+1) = gNamax.*(m(j,i).^3*h(j,i)).*(Vc(j,i)-ENa); 


        n(j,i+1) = 1-h(j,i);
        IK(j,i+1) = gKmax.*(n(j,i).^4).*(Vc(j,i)-Ek); 
        
      
        
        Ileak(j,i+1) = gleak.*(Vc(j,i)-Eleak); 
        Ids(j,i+1) = ((Vcb(j,i)-Vc(j,i))*(g/kap));
        

        

      nn=1.0;
   VA   = (((Vcb(j,i)-(Vc(j,i)))*g./(kap) + Ip(j,i) + nn*I(j,i) - gNamax*m(j,i)^3*h(j,i)*(Vc(j,i) - ENa) - gKmax*n(j,i).^4*(Vc(j,i) - Ek) - gleak*(Vc(j,i) - Eleak))/C)*dt;
   VB    = (((Vcb(j,i)-(Vc(j,i)+VA/2))*g./(kap) + Ip(j,i) + nn*I(j,i) - gNamax*m(j,i)^3*h(j,i)*((Vc(j,i)+ VA/2) - ENa) - gKmax*n(j,i).^4*((Vc(j,i) + VA/2) - Ek) - gleak*((Vc(j,i)+VA/2) - Eleak))/C)*dt;
   VC    = (((Vcb(j,i)-(Vc(j,i)+VB/2))*g./(kap) + Ip(j,i) + nn*I(j,i) - gNamax*m(j,i)^3*h(j,i)*((Vc(j,i)+ VB/2) - ENa) - gKmax*n(j,i).^4*((Vc(j,i) + VB/2) - Ek) - gleak*((Vc(j,i)+VB/2) - Eleak))/C)*dt;
   VD    = (((Vcb(j,i)-(Vc(j,i)+VC))*g./(kap)  + Ip(j,i) + nn*I(j,i) - gNamax*m(j,i)^3*h(j,i)*((Vc(j,i)+ VC) - ENa) - gKmax*n(j,i).^4*((Vc(j,i) + VC) - Ek) - gleak*((Vc(j,i)+VC) - Eleak))/C)*dt;
   Vc(j,i+1) = Vc(j,i) + (VA + 2*VB + 2*VC + VD)./6; 
        
        
   %-------------------------------------------------------DENDRITE    

        taumb(j,i+1) = (2*7.4/pi)*(26/(4*(Vcb(j,i)+45.67)^2 + 26^2));
         mbinf(j,i+1) = (-1/(1+exp((Vcb(j,i)+46.7)/5.7)) + 1);
  
         mb(j,i+1) = (mbinf(j,i)-((mbinf(j,i)-mb(j,i)).*(exp(-(dt)./(taumb(j,i))))));


 
 
        tauhb(j,i+1) = ((2*232/pi)*(43/(4*(Vcb(j,i)+60.0)^2 + 43^2)))*1.3;

          hbinf(j,i+1) = (1/(1+exp(-(Vcb(j,i)+55)/-3)));
        hb(j,i+1) = (hbinf(j,i)-((hbinf(j,i)-hb(j,i)).*(exp(-(dt)./(tauhb(j,i)))))); 
        INab(j,i+1) = gNamaxb.*(mb(j,i).^3*hb(j,i)).*(Vcb(j,i)-ENa); 
         

        taukb(j,i+1) = (0.40 + (2*70/pi)*(30/(4*(Vcb(j,i)+40)^2 + 30^2)));
        kbinf(j,i+1) = (1/(1+exp(-(Vcb(j,i)+12.5)/8.75)))^.25;
        kb(j,i+1) = (kbinf(j,i)-((kbinf(j,i)-kb(j,i)).*(exp(-(dt)./(taukb(j,i))))));

        
        
        IKv3b(j,i+1) = gKv3bmax.*(kb(j,i).^4).*(Vcb(j,i)-Ek); 

         Ileakb(j,i+1) = gleakb.*(Vcb(j,i)-Eleak); 


     
        bb=.0; %%0.5
   VAb    = (((Vc(j,i)-(Vcb(j,i)))*g./(1-kap) + bb*I(j,i)- gNamaxb*mb(j,i)^3*hb(j,i)*(Vcb(j,i) - ENa) - gKv3bmax*kb(j,i).^4*(Vcb(j,i) - Ek) - gKmaxb*nb(j,i).^4*(Vcb(j,i) - Ek) - gleakb*(Vcb(j,i) - Eleak))/Cb)*dt;
   VBb    = (((Vc(j,i)-(Vcb(j,i)+VAb/2))*g./(1-kap) + bb*I(j,i)- gNamaxb*mb(j,i)^3*hb(j,i)*((Vcb(j,i)+ VAb/2) - ENa) - gKv3bmax*kb(j,i).^4*((Vcb(j,i)+ VAb/2) - Ek) - gKmaxb*nb(j,i).^4*((Vcb(j,i) + VAb/2) - Ek) - gleakb*((Vcb(j,i)+VAb/2) - Eleak))/Cb)*dt;
   VCb    = (((Vc(j,i)-(Vcb(j,i)+VBb/2))*g./(1-kap)  + bb*I(j,i) - gNamaxb*mb(j,i)^3*hb(j,i)*((Vcb(j,i)+ VBb/2) - ENa) - gKv3bmax*kb(j,i).^4*((Vcb(j,i) + VBb/2) - Ek) - gKmaxb*nb(j,i).^4*((Vcb(j,i) + VBb/2) - Ek) - gleakb*((Vcb(j,i)+VBb/2) - Eleak))/Cb)*dt;
   VDb    = (((Vc(j,i)-(Vcb(j,i)+VCb))*g./(1-kap)  + bb*I(j,i) - gNamaxb*mb(j,i)^3*hb(j,i)*((Vcb(j,i)+ VCb) - ENa) - gKv3bmax*kb(j,i).^4*((Vcb(j,i) + VCb) - Ek) - gKmaxb*nb(j,i).^4*((Vcb(j,i) + VCb) - Ek) - gleakb*((Vcb(j,i)+VCb) - Eleak))/Cb)*dt;
   Vcb(j,i+1) = Vcb(j,i) + (VAb + 2*VBb + 2*VCb + VDb)./6; 
   
dVdt(j,i+1)=(Vc(j,i+1)-Vc(j,i))/dt;



  
    
    




end;




end;
figure(1);
plot (t,Vc);
figure(2);
plot(t,Vcb);
