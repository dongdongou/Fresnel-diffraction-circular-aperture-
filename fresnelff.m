clear;
clc;
step = 350;
lamda = 500e-6; %changed
k = 2*pi/lamda;
z = 12.5; %changed 
%ȷ��������
N = 500; %Բ����������
r = 0.25; %changed
I = zeros(N, N);
[m, n] = meshgrid(linspace(-N/step, N/step, N));
D = (m.^2+n.^2).^(1/2);
i = find(D <= r);
I(i) = 1;  %�հ뾶��Χ��͸��ϵ��Ϊ1
q = exp(j*k*(m.^2+n.^2)/2/z);
subplot(2,2,1); %Բ��ͼ��
imshow(I);
%imagesc(I) %������ͼ��
%colormap([0 0 0;1 1 1]) %�ڰ�����
% 
% I = I.*q;
L = 500;
M = 500; %ȡ��ͬ�������ھ�������
[x, y] = meshgrid(linspace(-L/step, L/step, M));
h = exp(j*k*z)*exp((j*k*(x.^2+y.^2))/(2*z))/(j*lamda*z); %������
%H = fftshift(fft2(h));
B = fftshift(fft2(I.*q));
G = h.*B; %
% U = fftshift(ifft2(G));
%Br = (abs(G)/max(abs(G))); %��һ��
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