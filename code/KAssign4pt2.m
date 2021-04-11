%Kwabena Gyasi Bawuah 
%101048814

close all 
clear

global G 
global C
global F
global Foff
global V2


%G = zeros(6, 6); 
%C = zeros(6, 6); 
F = zeros(7, 1);
Foff = zeros(7, 1);
Fsin = zeros(7,1);
Fg = zeros(7,1);


vo2 = zeros(1000, 1); 
W = zeros(1000, 1);
%part d
CC = zeros(1000,1);
GG = zeros(1000,1);


%Circuit Parameters 
%resistances and conductances
R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1; 
R0 = 1000; 
G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
G0 = 1/R0;
L = 0.2;
a = 100;
Cap = 0.25;
ts = 1000;
ts2 = 1.9898e4;
Vin = 1;
% Capacitance
Cn = 0.00001; 


vi = zeros(100, 1);
vo = zeros(100, 1);
v3 = zeros(100, 1);


V2 = zeros(7, ts);
V3 = zeros(7, ts);
V4 = zeros(7, ts2);

% G(1, 1) = 1;                                        
% G(2, 1) = G1; G(2, 2) = -(G1 + G2); G(2, 6) = -1;   
% G(3 ,3) = -G3; G(3, 6) = 1;                       
% G(4, 3) = -a*G3; G(4, 4) = 1;                        
% G(5, 5) = -(G4+G0); G(5, 4) = G4;   
% G(6, 2) = -1; G(6, 3) = 1;

G = [1 0 0 0 0 0 0;
           -G2 G1+G2 -1 0 0 0 0;
            0 1 0 -1 0 0 0;
            0 0 -1 G3 0 0 0;
            0 0 0 0 -a 1 0;
            0 0 0 G3 -1 0 0;
            0 0 0 0 0 -G4 G4+G0];


% C(2, 1) = Cap; C(2, 2) = -Cap;
% C(6, 6) = L;

C = [0 0 0 0 0 0 0;
           -Cap Cap 0 0 0 0 0;
            0 0 -L 0 0 0 0;
            0 0 0 -Cn 0 0 0;
            0 0 0 0 0 0 0;
            0 0 0 -Cn 0 0 0;
            0 0 0 0 0 0 0;];

      

v = 0;
% [V1 Iin V2 I3 V4 Icc Vo]

%DC Case
%capacitive/inductive effect
s = 0;


F(1,1) = Vin;
 
Foff(1,1) = Vin-Vin;


%2D:ii:A 0.03s time step
V1 = zeros(7, ts);
Vstart = zeros(7, 1);
dt=1e-3;
dt2 = 1.9898e-4;


for i = 1:ts
    
    F(1,1) = exp(-1/2*((i/ts-0.06)/(0.03))^2);
    F(4,1) = 0.001*randn();
    F(7,1) = 0.001*randn();
    
    if i == 1
        V1(:,i) = (C./dt+G)\(F+C*Vstart/dt);
    else
        V1(:,i) = (C./dt+G)\(F+C*Vpast/dt);
    end
    Vpast = V1(:, i);
end

figure(1)
subplot(2,1,1)
plot(1:ts, V1(7,:))
hold on
plot(1:ts, V1(1,:))
title('Vout with Noise')
xlabel('Time(ms)')
ylabel('V(v)')
grid on

% Frequency domain
f = (-ts/2:ts/2-1);              

fV1 = fft(V1.');
fsV1 = fftshift(fV1);
%simulation
figure(1)
subplot(2,1,2)
plot(f, abs(fsV1(:, 1)))
hold on
plot(f, abs(fsV1(:, 7)))
title('Fourier-Transform of Vout')
xlabel('w(1/ms)')
ylabel('V(v)')
grid on

for i = 1:ts
    
    F(1,1) = exp(-1/2*((i/ts-0.06)/(0.03))^2);
    F(4,1) = 0.001*randn();
    F(7,1) = 0.001*randn();
    if i == 1
        V3(:,i) = (C./dt+G)\(F+C*Vstart/dt);
    else
        V3(:,i) = (C./dt+G)\(F+C*Vpast/dt);
    end
    Vpast = V3(:, i);
        
end
figure(2)
plot(1:ts, V3(7,:))
hold on
plot(1:ts, V3(1,:))
title('Vin and Vout with initial Cn')
xlabel('Time(ms)')
ylabel('V(v)')
grid on

C(4,4) = -1e-13;
C(6,4) = -1e-13;

for i = 1:ts
    
    F(1,1) = exp(-1/2*((i/ts-0.06)/(0.03))^2);
    F(4,1) = 0.001*randn();
    F(7,1) = 0.001*randn();
    if i == 1
        V2(:,i) = (C./dt+G)\(F+C*Vstart/dt);
    else
        V2(:,i) = (C./dt+G)\(F+C*Vpast/dt);
    end
    Vpast = V2(:, i);
        
end
figure(3)
subplot(2,1,1)
plot(1:ts, V2(7,:))
hold on
plot(1:ts, V2(1,:))
title('Vin and Vout with small Cn')
xlabel('Time(ms)')
ylabel('V(v)')
grid on

C(4,4) = -5.3e-5;
C(6,4) = -5.3e-12;

for i = 1:ts
    
    F(1,1) = exp(-1/2*((i/ts-0.06)/(0.03))^2);
    F(4,1) = 0.001*randn();
    F(7,1) = 0.001*randn();
    if i == 1
        V3(:,i) = (C./dt+G)\(F+C*Vstart/dt);
    else
        V3(:,i) = (C./dt+G)\(F+C*Vpast/dt);
    end
    Vpast = V3(:, i);
        
end
figure(3)
subplot(2,1,2)
plot(1:ts, V3(7,:))
hold on
plot(1:ts, V3(1,:))
title('Vin and Vout with large Cn')
xlabel('Time(ms)')
ylabel('V(v)')
grid on

%%%%%%%%%%%%%%

C(4,4) = -Cn;
C(6,4) = -Cn;

for i = 1:ts2
    
    Fg(1,1) = exp(-1/2*((i/ts2-0.06)/(0.03))^2);
    Fg(4,1) = 0.001*randn();
    Fg(7,1) = 0.001*randn();
    if i == 1
        V4(:,i) = (C./dt2+G)\(Fg+C*Vstart/dt2);
    else
        V4(:,i) = (C./dt2+G)\(Fg+C*Vpast/dt2);
    end
    Vpast = V4(:, i);
        
end
figure(4)
plot(1:ts2, V4(7,:))
hold on
plot(1:ts2, V4(1,:))
title('Vin vs Vout with change in time steps')
xlabel('Time(ms)')
ylabel('V(v)')
grid on