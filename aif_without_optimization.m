clc
clear all
close all







%% 场群偏航设置
matrix=zeros(1,159);


%% 场群功率优化
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

swi=SmartWindInterface_aif(turbine_diameter_vector,turbine_hub_height_vector,matrix);
rng("default")
% swi.windfield.wind_direction=270;
swi.windfield.wake.velocity_model='Gauss';
swi.windfield.wake.deflection_model='Gauss';
swi.windfield.enable_wfr='yes';
swi.windfield.resolution=[10 10 20];

swi.calculate_wake();
p1=swi.get_farm_power();
y1=swi.get_aif_angles();
figure(1)


swi.show_horplane(140);

a=[];
b=[];
for i=1:size(swi.layout_x)
 if swi.windfield.turbinechart.turbines{i,1}.hub_height==108
    a=[a,i];
 else
    b=[b,i];

 end

end

p_68=[];
p_83=[];
for i=1:size(a,2)


p_68=[p_68,swi.windfield.turbinechart.turbines{a(i),1}.power];




end
sum_a=sum(p_68);

for i=1:size(b,2)


p_83=[p_83,swi.windfield.turbinechart.turbines{b(i),1}.power];




end
sum_b=sum(p_83);





























