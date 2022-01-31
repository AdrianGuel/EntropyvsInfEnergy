%Entropy production self-assembly line formation of particles
%Author Adrian Guel December 2021

clear all;
close all;

N=3;
a=0;
b=5;
rng(1);
kp=(a + (b-a)*rand)*ones(N,1);
rng(4);
m=a + (b-a).*rand(N,1);


A=zeros(2*N,2*N);
S0=eye(2*N,2*N).*[kp/5; kd/5];
rng(5);
X0=a + (b-a).*rand(2*N,1);