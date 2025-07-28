clc
clear all
close all

tic

swi=SmartWindInterface();
swi.windfield.wake.velocity_model='Gauss';
swi.windfield.wake.deflection_model='Gauss';
swi.windfield.enable_wfr='yes';
swi.windfield.resolution=[500 300 100];
swi.calculate_wake();
figure
swi.show_horplane(90);
figure
swi.show_verplane(10);
figure
swi.show_crossplane(630);

toc



























