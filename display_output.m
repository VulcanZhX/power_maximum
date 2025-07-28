clc
clear all
close all

load swi.mat
load p1.mat
load ans.mat
load y1.mat

swi.calculate_wake();


% swi=SmartWindInterface();
% swi.windfield.wake.velocity_model='Gauss';
% swi.windfield.wake.deflection_model='Gauss';
% swi.windfield.enable_wfr='yes';
% swi.windfield.resolution=[250 150 50];
% swi.calculate_wake();
% figure(1)
% swi.show_horplane(140);



figure


swi.show_horplane(140);
























































