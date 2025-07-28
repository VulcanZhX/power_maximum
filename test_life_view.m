clc
clear all
close all





load('size_50_50_50.mat')


obj=swi;

for i=1:size(obj.windfield.turbinechart.layout_x)





     life_single_turbine(i,1)=obj.windfield.turbinechart.turbines{i,1}.left_life;
     turbulence_single_turbine(i,1)=obj.windfield.turbinechart.turbines{i,1}.turbulence;





end

a=mean(life_single_turbine);
b=mean(turbulence_single_turbine); 

a_max=max(life_single_turbine);
a_min=min(life_single_turbine);
aa=sqrt(mean((life_single_turbine-mean(life_single_turbine)).^2));































