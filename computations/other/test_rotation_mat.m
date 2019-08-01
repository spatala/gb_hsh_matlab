%#ok<*FLUDLR>

% rng(142857);

w  = 8 * pi * (rand() - 0.5);
th = pi * rand();
ph = 2 * pi * rand();

j = 3;

% Explicit form for j = 1

% ref = [(cos(w / 2) - 1i * sin(w / 2) * cos(th))^2, ...
%        -1i * sqrt(2) * sin(w / 2) * sin(th) * exp(-1i * ph) * (cos(w / 2) - 1i * sin(w / 2) * cos(th)), ...
%        -(sin(w / 2) * sin(th) * exp(-1i * ph))^2;
%        -1i * sqrt(2) * sin(w / 2) * sin(th) * exp(1i * ph) * (cos(w / 2) - 1i * sin(w / 2) * cos(th)), ...
%        1 - 2 * sin(w / 2)^2 * sin(th)^2, ...
%        -1i * sqrt(2) * sin(w / 2) * sin(th) * exp(-1i * ph) * (cos(w / 2) + 1i * sin(w / 2) * cos(th));
%        -(sin(w / 2) * sin(th) * exp(1i * ph))^2, ...
%        -1i * sqrt(2) * sin(w / 2) * sin(th) * exp(1i * ph) * (cos(w / 2) + 1i * sin(w / 2) * cos(th)), ...
%        (cos(w / 2) + 1i * sin(w / 2) * cos(th))^2];
% ref = flipud(fliplr(ref));

% U = rotation_mat(j, [w-2*pi, th, ph]);
U = rotation_mat(j, w-2*pi, th, ph);

R = rotation([cos(ph) * sin(th), sin(ph) * sin(th), cos(th)], w, 2 * j);
R = flipud(fliplr(R));

disp(['Max U R difference: ', num2str(max(abs(U(:) - R(:))))]);

tic;
for a = 1:1000
%     U = rotation_mat(j, [w, th, ph]);
    U = rotation_mat(j, w, th, ph);
end
toc;

ax = [cos(ph) * sin(th), sin(ph) * sin(th), cos(th)];
tic;
for a = 1:1000
    R = rotation(ax, w, 2 * j);
end
toc;