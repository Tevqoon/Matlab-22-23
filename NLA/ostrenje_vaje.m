% to poskrbi, da je X vedno ista random matrika

nakljucno = rng('default');
rng(2);

% nalozimo originalno sliko
% X = imread('butterfly.jpg');
X = imread('pumpkins.tif');
% X = imread("mimi.jpeg");

% crnobela tehnika
% X = rgb2gray(X);
X = im2double(X);
[n, m] = size(X);
%fig = figure;
fig = 1;
subplot(2,3,1);
imagesc(X), axis image, colormap(gray);







% parameter s = velikost Gaussove zameglitve

s = 2;

[PSF, center] = psfGauss([n,m],s);

[Ar, Ac] = kronDecomp(PSF, center);





% zameglimo slike z leve in desne

B = Ac * X * Ar';

subplot(2,3,2);

imagesc(B), axis image, colormap(gray);





% parameter e = velikost suma

e = 0.1;

E = e * randn(n,m);

Bhat = B + E ; % Bhat je zamegljena slika, kateri smo dodali tudi sum E

subplot(2,3,3);

imagesc(Bhat), axis image, colormap(gray);





% SVD za Ar in Ac

[Ur, Sr, Vr] = svd(Ar);

[Uc, Sc, Vc] = svd(Ac);

W = Sc\Uc'* Bhat * Ur/Sr;
[wx, wy] = size(W);

% naiven nacin brez regularizacije

naivniX = Vc * (Sc\Uc'* Bhat * Ur/Sr) * Vr';

subplot(2,3,4);

imagesc(naivniX), axis image, colormap(gray);





% odrezani SVD, glej 0.2.3

odrez = 30;

sigma = diag(Sc) * diag(Sr)';

%.................

t = 0.35;
F_odr = sigma >= t;
M_odr = F_odr .* W;

X_odrezanSVD = Vc * M_odr * Vr';



subplot(2,3,5);

imagesc(X_odrezanSVD), axis image, colormap(gray);


% SVD Tihonov, glej 0.2.3


alpha = 100;

Fi = (sigma .^ 2) ./ (sigma .^ 2 + alpha^2);

M_Tih = Fi .* W;

X_Tihonov = Vc * M_Tih * Vr';

figure(fig)

subplot(2,3,6);

imagesc(X_Tihonov), axis image, colormap(gray);