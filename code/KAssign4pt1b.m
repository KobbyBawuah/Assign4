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

vi = zeros(100, 1);
vo = zeros(100, 1);
v3 = zeros(100, 1);


V2 = zeros(7, ts);
V3 = zeros(7, ts);

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
            0 0 0 0 0 0 0;
            0 0 0 0 0 0 0;
            0 0 0 0 0 0 0;
            0 0 0 0 0 0 0;];


v = 0;
% [V1 Iin V2 I3 V4 Icc Vo]

%DC Case
%capacitive/inductive effect
s = 0;

Vin = 1;

F(1,1) = Vin;
 
Foff(1,1) = Vin-Vin;


%2D:ii:A 0.03s time step
V1 = zeros(7, ts);
Vstart = zeros(7, 1);
dt=1e-3;

for i = 1:ts
    if i < 30
        V1(:,i) = (C./dt+G)\(Foff+C*Vstart/dt);
    elseif i == 30
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
title('step that transitions 0-1 at 0.03s')
xlabel('Time(ms)')
ylabel('V(v)')
grid on

for i = 1:ts

    Vsin = sin(2*3.14*(1/0.03)*i/ts);
    Fsin(1,1) = Vsin;
    if i == 1
        V2(:,i) = (C./dt+G)\(Fsin+C*Vstart/dt);
    else
        V2(:,i) = (C./dt+G)\(Fsin+C*Vpast/dt);
    end
    Vpast = V2(:, i);   
end

figure(1)
subplot(2,1,2)
plot(1:ts, V2(7,:))
hold on
plot(1:ts, V2(1,:))
title('sin function Vin and Vo')
xlabel('Time(ms)')
ylabel('V(v)')
grid on


for k = 1:ts

    Vgauss = exp(-1/2*((k/ts-0.06)/(0.03))^2);
    Fg(1,1) = Vgauss;
    if k == 1
        V3(:,k) = (C./dt+G)\(Fg+C*Vstart/dt);
    else
        V3(:,k) = (C./dt+G)\(Fg+C*Vpast/dt);
    end
    Vpast = V3(:, k);
        
end

figure(2)
plot(0:ts-1, V3(7,:))
hold on
plot(0:ts-1, V3(1,:))
title('Vin and Vout: gaussian function')
xlabel('T(ms)')
ylabel('V(v)')
grid on

%Plot of Vin, Vo withinput signal in f-domain
f = (-ts/2:ts/2-1); 
fV1in = fft(V1(1, :));
fV1out = fft(V1(7, :));
fsV1in = fftshift(fV1in);
fsV1out = fftshift(fV1out);

figure(3)
plot(f, abs(fsV1in))
hold on
plot(f, abs(fsV1out))
xlim([-150,150]);
title('Vin and Vout: frequency domain (step function)')
xlabel('w(1/ms)')
ylabel('V(v)')
grid on

%sine function 
fV2 = fft(V2.');
fsV2 = fftshift(fV2);
figure(4)
plot(f, abs(fsV2(:, 1)))
hold on
plot(f, abs(fsV2(:, 7)))
xlim([-150,150]);
title('Vin and Vout: frequency domain (sin function)')
xlabel('w')
ylabel('V(v)')

%guass function
fV3 = fft(V3.');
fsV3 = fftshift(fV3);
figure(5)
plot(f, abs(fsV3(:, 1)))
hold on
plot(f, abs(fsV3(:, 7)))
xlim([-150,150]);
title('Vin and Vout: frequency domain (gauss function)')
xlabel('w')
ylabel('V(v)')


