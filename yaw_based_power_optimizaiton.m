clc
clear all
close all







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEST_10%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sqz_1=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','WindField','Range','B10:B10'));              %青州1风机个数
sqz_2=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','WindField','Range','B11:B11'));              %青州2风机个数
sqz_3=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','WindField','Range','B12:B12'));              %青州3风机个数
turbine_diameter_vector=[cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B3:B3')),...
                cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','E3:E3')),...
                cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','H3:H3'))];
turbine_hub_height_vector=[cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B4:B4')),...
                cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','E4:E4')),...
                cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','H4:H4'))];
sqz_12=sqz_1+sqz_2;

swi=SmartWindInterface(sqz_12,turbine_diameter_vector,turbine_hub_height_vector);
rng("default")
% swi.windfield.wind_direction=270;
swi.windfield.wake.velocity_model='Gauss';
swi.windfield.wake.deflection_model='Gauss';
swi.windfield.enable_wfr='yes';
swi.windfield.resolution=[100 100 20];

swi.calculate_wake();
p1=swi.get_farm_power();
y1=swi.get_yaw_angles();
figure(1)


swi.show_horplane(140);

tic
swi.yaw_optimization_gb();
toc
p2=swi.get_farm_power();
y2=swi.get_yaw_angles();
swi.calculate_wake();
figure(2)


swi.show_horplane(140);



optimized=(p2-p1)/p1;
% swi.reset_farm_keep_layout();
% swi.yaw_optimization_ga();
% p2=swi.get_farm_power();
% y2=swi.get_yaw_angles();
% swi.reset_farm_keep_layout();
% tic
% swi.yaw_optimization_sq();
% toc
% p3=swi.get_farm_power();
% y3=swi.get_yaw_angles();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEST_11%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% swi=SmartWindInterface();
% rng("default")
% layout=generate_random_layout_array(8,1400,1400,500);
% swi.set_layout(layout);
% swi.ya_options=6;
% swi.yaw_optimization_miga();
% p4=swi.get_farm_power();
% y4=swi.get_yaw_angles();
% swi.reset_farm_keep_layout();
% swi.ya_options=6;
% swi.yaw_optimization_sq();
% p5=swi.get_farm_power();
% y5=swi.get_yaw_angles();
% 
% swi.calculate_wake();
% swi.show_horplane(90);
% swi.windfield.wind_direction
% swi.windfield.turbinechart.turbines





























