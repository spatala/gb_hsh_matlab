% rng(142857);

w1  = pi * rand();
th1 = pi * rand();
ph1 = 2. * pi * rand();

CK = angles_to_CK(w1, th1, ph1);
[w2, th2, ph2] = CK_to_angles(CK);

disp(['Max difference: ', num2str([w1, th1, ph1] - [w2, th2, ph2])]);
