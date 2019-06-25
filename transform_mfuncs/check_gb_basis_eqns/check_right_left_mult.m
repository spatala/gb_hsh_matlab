clear all; clc;
nsz = 3;

A1 = rand(nsz,nsz); 
A1v = reshape(A1', [1,nsz*nsz]); 

Ar = rand(nsz,nsz); 
Al = rand(nsz,nsz);

% save('A_mats','A1', 'Ar', 'Al');
% s1 = load('A_mats.mat');
% A1 = s1.A1; Ar = s1.Ar; Al = s1.Al;

% Amult = Al*A1;
% Amult_v1 = reshape(Amult', [1,nsz*nsz]); 
% Amult_v2 = A1v*kron(Al',eye(nsz,nsz));

Amult = A1*Ar;
Amult_v1 = reshape(Amult', [1,nsz*nsz]); 
Amult_v2 = A1v*kron(eye(nsz,nsz),Ar);

% Amult = Al*A1*Ar;
% Amult_v1 = reshape(Amult', [1,nsz*nsz]); 
% Amult_v2 = A1v*kron(Al',Ar);

norm(Amult_v1 - Amult_v2)