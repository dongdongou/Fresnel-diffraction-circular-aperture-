clear;
clc;
step = 350;
lamda = 500e-6; %changed
k = 2*pi/lamda;
z = 12.5; %changed 
%确定衍射屏
N = 500; %圆屏采样点数
r = 0.25; %changed
I = zeros(N, N);
[m, n] = meshgrid(linspace(-N/step, N/step, N));
D = (m.^2+n.^2).^(1/2);
i = find(D <= r);
I(i) = 1;  %空半径范围内透射系数为1
q = exp(j*k*(m.^2+n.^2)/2/z);
subplot(2,2,1); %圆孔图像
imshow(I);
%imagesc(I) %衍射屏图像
%colormap([0 0 0;1 1 1]) %黑白区分
% 
% I = I.*q;
L = 500;
M = 500; %取相同点数用于矩阵运算
[x, y] = meshgrid(linspace(-L/step, L/step, M));
h = exp(j*k*z)*exp((j*k*(x.^2+y.^2))/(2*z))/(j*lamda*z); %接收屏
%H = fftshift(fft2(h));
B = fftshift(fft2(I.*q));
G = h.*B; %
% U = fftshift(ifft2(G));
%Br = (abs(G)/max(abs(G))); %归一化
C = abs(G);
subplot(2,2,2);imagesc(C);
axis image;
colormap(hot);
% %figure;
subplot(2,2,3);mesh(x,y,abs(G));
subplot(2,2,4);
axis image;
d = C(251,:);
d = d/max(d);
plot(d);