function dSt=rasisM1(t,St,p)
%Differential equations for the model without vaccination

global B wm wi1 wi2 u um d1 d2 rr1 rr2 ri2 ri3 v1 v2 v3 reintro wA;  

b=p(1);
b1=p(2);
phi1=p(3);
al=p(4);
ar1=1; %p(5);
ar2=1; %p(6);

%c1=100*ones(al); %Homogeneous mixing
c2=100*[ar1*ones(al,12) ar2*ones(al,12) ones(al,al-24)]; %Age-specific risk of acquisition
beta=(b/100)*c2; %/sum(N) type of mixing

lamda=(1+b1*cos(2*pi*(t-phi1)/52.18))*beta*(St(2*al+1:3*al)+ri2*St(5*al+1:6*al)+ri3*St(8*al+1:9*al))/sum(St);

for i=1:al
    dSt(i,1)=log(1+B(round(t),i))/52*sum(St(:))*(1-v1(round(t),1))  - (u(i)+um(i))*St(i) - wm*St(i); %dM/dt  
    dSt(i+al,1)=wm*St(i+14*al) - St(i+al)*(sum(lamda(i,:))+reintro) - (u(i)+um(i))*St(i+al); %dS0/dt 
    dSt(i+2*al,1)=St(i+al)*(sum(lamda(i,:))+reintro) - (d1+u(i)+um(i))*St(i+2*al); %dI1/dt
    dSt(i+3*al,1)=d1*St(i+2*al) - (wi1+u(i)+um(i))*St(i+3*al) + v1(round(t),1)*log(1+B(round(t),i))/52*sum(St(:)); %dR1/dt
    dSt(i+4*al,1)=wi1*St(i+3*al) - St(i+4*al)*rr1*(sum(lamda(i,:))+reintro) - (u(i)+um(i))*St(i+4*al) + wA*St(i+7*al); %dS1/dt 
    dSt(i+5*al,1)=St(i+4*al)*rr1*(sum(lamda(i,:))+reintro) - (d2+u(i)+um(i))*St(i+5*al); %dI2/dt
    dSt(i+6*al,1)=d2*St(i+5*al) - (wi1+u(i)+um(i))*St(i+6*al); %dR2/dt
    dSt(i+7*al,1)=wi1*St(i+6*al) + wi2*St(i+9*al) - St(i+7*al)*rr2*(sum(lamda(i,:))+reintro) - (u(i)+um(i)+wA)*St(i+7*al); %dSR/dt
    dSt(i+8*al,1)=St(i+7*al)*rr2*(sum(lamda(i,:))+reintro) - (d2+u(i)+um(i))*St(i+8*al); %dAI/dt
    dSt(i+9*al,1)=d2*St(i+8*al) - (wi2+u(i)+um(i))*St(i+9*al); %dR/dt
    
    dSt(i+10*al,1)=wm*St(i) - (u(i)+um(i))*St(i+10*al) - wm*St(i+10*al); %dM2/dt
    dSt(i+11*al,1)=wm*St(i+10*al) - (u(i)+um(i))*St(i+11*al) - wm*St(i+11*al); %dM3/dt
    dSt(i+12*al,1)=wm*St(i+11*al) - (u(i)+um(i))*St(i+12*al) - wm*St(i+12*al); %dM4/dt
    dSt(i+13*al,1)=wm*St(i+12*al) - (u(i)+um(i))*St(i+13*al) - wm*St(i+13*al); %dM5/dt
    dSt(i+14*al,1)=wm*St(i+13*al) - (u(i)+um(i))*St(i+14*al) - wm*St(i+14*al); %dM6/dt
    
    if i>1 %aging of individuals from one age group to the next and vaccination
        dSt(i,1)=dSt(i,1) + (1-v1(round(t),i))*u(i-1)*St(i-1);
        dSt(i+al,1)=dSt(i+al,1) + (1-v1(round(t),i))*u(i-1)*St(i+al-1);
        dSt(i+2*al,1)=dSt(i+2*al,1) + u(i-1)*St(i+2*al-1);
        dSt(i+3*al,1)=dSt(i+3*al,1) + (1-v2(round(t),i))*u(i-1)*St(i+3*al-1) + v1(round(t),i)*u(i-1)*(St(i-1)+St(i+al-1)+St(i+10*al-1)+St(i+11*al-1)+St(i+12*al-1)+St(i+13*al-1)+St(i+14*al-1));
        dSt(i+4*al,1)=dSt(i+4*al,1) + (1-v2(round(t),i))*u(i-1)*St(i+4*al-1); 
        dSt(i+5*al,1)=dSt(i+5*al,1) + u(i-1)*St(i+5*al-1);
        dSt(i+6*al,1)=dSt(i+6*al,1) + (1-v3(round(t),i))*u(i-1)*St(i+6*al-1) + v2(round(t),i)*u(i-1)*(St(i+3*al-1)+St(i+4*al-1));
        dSt(i+7*al,1)=dSt(i+7*al,1) + (1-v3(round(t),i))*u(i-1)*St(i+7*al-1);
        dSt(i+8*al,1)=dSt(i+8*al,1) + u(i-1)*St(i+8*al-1);
        dSt(i+9*al,1)=dSt(i+9*al,1) + u(i-1)*St(i+9*al-1) + v3(round(t),i)*u(i-1)*(St(i+6*al-1)+St(i+7*al-1));
                
        dSt(i+10*al,1)=dSt(i+10*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+10*al-1);
        dSt(i+11*al,1)=dSt(i+11*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+11*al-1);
        dSt(i+12*al,1)=dSt(i+12*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+12*al-1);
        dSt(i+13*al,1)=dSt(i+13*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+13*al-1);
        dSt(i+14*al,1)=dSt(i+14*al,1) + (1-v1(round(t),i))*u(i-1)*St(i+14*al-1);    
        
    end
end
