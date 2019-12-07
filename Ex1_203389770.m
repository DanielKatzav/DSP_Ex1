clear all
close all
clc

% ----------------- Question 1 ----------------- %

w_m = 3*pi;
dt = 1/100;
t = 0.2:dt:3;

x_t = 8./(w_m*t.^2) .* ((sin((1/2)*w_m*t)).^3).*cos(2*w_m*t);

figure(1)
plot(t,x_t)
title('Signal x(t) - Time domain')
xlabel('Time [sec]')
ylabel('x(t) [V]')

dw = 0.1;
w = -17*pi:dw:17*pi;

w1 = pi.*(w-w_m/2).*sign(w-w_m/2);
w2 = pi.*(w+w_m/2).*sign(w+w_m/2);
w3 = pi.*(w-(3*w_m)/2).*sign(w-(3*w_m)/2);
w4 = pi.*(w+(3*w_m)/2).*sign(w+(3*w_m)/2);

X1_w = (1/(2i)).*(3.*(pi.*(-w_m/2).*sign(-w_m/2)-pi.*(+w_m/2).*sign(+w_m/2)) - (pi.*(-(3*w_m)/2).*sign(-(3*w_m)/2)-pi.*(+(3*w_m)/2).*sign(+(3*w_m)/2)));

X2_w = -2i*(3*(WindowFunction(w,(2*pi/w_m), (w_m/2)) - WindowFunction(w,(2*pi/w_m),-w_m/2)) - ( WindowFunction(w,(2*pi/w_m), (3*w_m/2))- WindowFunction(w,(2*pi/w_m), (-3*w_m/2))));


Xw = X1_w + X2_w;

figure(2)
plot(w,abs(Xw))
title('Signal |X(j\omega)| - Scaled Frequency domain')
xlabel('\omega [rad/sec]')
ylabel('X(j\omega) [V]')
xticklabels({'-17\pi','-10\pi','-5\pi','0','5\pi','10\pi','-17\pi'})
xlim([-17*pi,17*pi])

w_max = 4.6441*pi;
w_s = 8*w_max;
Ts = ((2*pi)/w_s);
dn = Ts/length(x_t);
n = 0:dn:Ts-dn;
x_zoh = ZOH_me(t,x_t,Ts);

figure(3)
plot(t,x_t,'k','LineWidth',2)
hold on
plot(t,x_zoh)
legend([{'x(t)'};{'x_Z_O_H(t)'}])
title('Question 1.c - x(t) and x_Z_O_H in Time [sec]')
xlabel('Time[sec]')
ylabel('Signal x(t) && x_Z_O_H(t) in [V] - 1.c')

X_ZOH = fft(x_zoh);
figure(4)
plot(abs(X_ZOH))





function x_ZOH = ZOH_me(t, x_t, Ts)
    %ZOH Sampling of x(t) with intervals n*Ts 
    x_ZOH = zeros(1,length(x_t));

    for n = 1:(length(t)-1)*Ts
        x_ZOH(1+10*(n-1):10*n) = x_ZOH(1+10*(n-1):10*n) + x_t(1+10*(n-1));
    end
end
function res = WindowFunction(length_vector, region_length, region_start)
    dw = 0.1;

    starting_point = round(region_start/dw +round(length(length_vector)/2));
    res = zeros(length(length_vector),1);
    res(starting_point:starting_point+region_length,1) = 1;
   
end