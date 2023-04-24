format long

disp("Prva naloga")
a1 = 20316 / 1e4;
p1 = [1, 6, -5, -30, 4, 24];

disp(horner_with_float(p1, a1, 5, 10))
fl1 = @(x) float(x, 5, 10);
p1x = fl1((a1 - 2) * fl1((a1 - 1) * fl1((a1 + 1) * fl1((a1 + 2) * (a1 + 6)))));
disp(p1x)

disp("Druga naloga")
b2 = 50/199;
vir = [b2, 0.5 * b2];
smer0 = [1, -1]';
[t1, s1] = odboj_druga_test(vir, smer0);
disp(t1(1) + t1(2))
[t2, s2] = odboj_druga_test(t1, s1);
[t3, s3] = odboj_druga_test(t2, s2);
[t4, s4] = odboj_druga_test(t3, s3);

pot = norm(vir - t1) + norm(t1 - t2) + norm(t3 - t2) + norm(t4 - t3);

disp(pot)

disp("Tretja naloga")

c3 = 55 / 501;
f3 = @(x) sin(exp(-9 * x^2)) - c3;
[y4, korr31] = regula_falsi(f3, 0, 1, 0, 4);
disp(y4)
[yn, korr32] = regula_falsi(f3, 0, 1, 1e-3, Inf);
disp(yn)

disp("Četrta naloga")

d4 = -132/101;
% A = K_cetrta([1 2 3 4 5 6 7 8 9 d4]);
A = [1 2 3 4 5 6 7 8 9 d4; 
    2 3 4 5 6 7 8 9 d4 1; 
    3 4 5 6 7 8 9 d4 1 2; 
    4 5 6 7 8 9 d4 1 2 3;
    5 6 7 8 9 d4 1 2 3 4;
    6 7 8 9 d4 1 2 3 4 5;
    7 8 9 d4 1 2 3 4 5 6;
    8 9 d4 1 2 3 4 5 6 7;
    9 d4 1 2 3 4 5 6 7 8;
    d4 1 2 3 4 5 6 7 8 9];
disp(norm(A, "fro") / norm(A, 2))

[L, U, P] = lu(A);
Pb = P * (d4 * ones(10, 1));
y = U \ Pb;

disp(norm(y, "inf"))