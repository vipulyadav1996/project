function [ssv,w_n] = FRF_Cb(w,F0,k2,c2,nod)

t_time=800;
t(1)=0;
x0=0;y0=0;u0=0;w0=0;

m1=50;m_abs=5;
k1=11250; 
%k2=2.347; 
k3=14238; %k3=non-linearity in primary system(1000N/m^3)
c1=30;
k4=0;
%k4=7.32;           %4.53;
%N = m_abs*9.81;
%del = m_abs*9.81*2  ;  %friction force
%0.11;
w_f=w;            
            
w_n=(k1/m1)^0.5;

x(1) = x0 ; y(1) = y0 ; u(1) =  u0 ; w(1) = w0 ; %initial condition x-mass,y-absorber
n=30000;


f1= @(x, y, u, w, t) (F0*sin(w_f*t)-k1*x-k3*x^3-k2*(x-y)-c1*(u)-c2*(u-w)-k4*(x-y)^3)/m1;
f2 = @(x, y, u, w) (-c2*(w-u)-k2*(y-x)-k4*(y-x)^3)/m_abs;


h = t_time / n;
 
       
for i = 1 : n
            %update time
            t(i+1)=t(i)+h;
               
            %RK method 
                kx1 =  u(i);
                ky1 = w(i);
                ku1 = f1(x(i), y(i), u(i), w(i),t(i));
                kw1 = f2(x(i), y(i), u(i), w(i));
               
                kx2 = (u(i) + h*ku1/2);
                ky2 = (w(i) + h*kw1/2);
                ku2 = f1(x(i)+h*kx1/2, y(i)+h*ky1/2, u(i)+h*ku1/2, w(i)+h*kw1/2, t(i)+h/2);
                kw2 = f2(x(i)+h*kx1/2, y(i)+h*ky1/2, u(i)+h*ku1/2, w(i)+h*kw1/2);
                
                kx3 = (u(i) + h*ku2/2);
                ky3 = (w(i) + h*kw2/2);
                ku3 = f1(x(i)+h*kx2/2, y(i)+h*ky2/2, u(i)+h*ku2/2, w(i)+h*kw2/2, t(i)+h/2);
                kw3 = f2(x(i)+h*kx2/2, y(i)+h*ky2/2, u(i)+h*ku2/2, w(i)+h*kw2/2);
               
                kx4 = (u(i) + h*ku3);
                ky4 = (w(i) + h*kw3);
                ku4 = f1(x(i)+h*kx3, y(i)+h*ky3, u(i)+h*ku3, w(i)+h*kw3, t(i)+h);
                kw4 = f2(x(i)+h*kx3, y(i)+h*ky3, u(i)+h*ku3, w(i)+h*kw3);
               
                
               
                x(i+1) = x(i) + h*(kx1 + 2*kx2 + 2*kx3 + kx4) / 6;
                y(i+1) = y(i) + h*(ky1 + 2*ky2 + 2*ky3 + ky4) / 6;
                u(i+1) = u(i) + h*(ku1 + 2*ku2 + 2*ku3 + ku4) / 6;
                w(i+1) = w(i) + h*(kw1 + 2*kw2 + 2*kw3 + kw4) / 6;
                
               
end

peak_local=findpeaks(x);
[r c] = size(peak_local);
if nod == 'y'
    plot(t,x)
end
ssv = peak_local(c-2)/(F0/k1);
end


