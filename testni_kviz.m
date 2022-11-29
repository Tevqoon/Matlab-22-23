format long

disp("Prva naloga")
p1 = [1 -8 23 -28 12];
a1 = 4.21;

disp(horner_with_float(p1, a1, 52, 2))
disp(horner_with_float(p1,a1, 5, 10))
disp(horner_with_float(p1,a1,6,10))

disp("Druga naloga")
vir = [3.5, 0.5];
smer0 = [1 2]';
[t1, s1, i1] = odboj_druga(vir, smer0);
[t2, s2, i2] = odboj_druga(t1, s1);
[t3, s3, i3] = odboj_druga(t2, s2);
[t4, s4, i4] = odboj_druga(t3, s3);

disp([i1 i2 i3 i4])
disp(t1(1)) 
disp(t2(2)) 
disp(norm(t4 - vir))

disp("Tretja naloga")

f3 = @(x) 5 * cos(x - exp(x)) - x;
[y21, korr21] = regula_falsi(f3, 2.75, 3, 0, 3);
disp(y21)
[y22, korr22] = regula_falsi(f3, 2, 2.5, 1e-10, Inf);
disp(korr22)
[y23, korr23] = regula_falsi(f3, 0, 3, 0, Inf);
disp(y23)

disp("Četrta naloga")

a4 = primes(100);
K = K_cetrta(a4);
l4 = length(a4);

disp(norm(K,1) + norm(K, 2) + norm(K, "inf") + norm(K, "fro"))
[Lb, Ub] = lubp(K);
disp(norm(Lb\ones(l4, 1), 2))
disp(norm(K\ones(l4, 1), 1))

[Ld, Ud] = lu(K);
disp(norm(Ub(:), "inf") / norm(Ud(:), "inf"))