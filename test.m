clc;
clear all;
close all;


load('swi.mat')


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



