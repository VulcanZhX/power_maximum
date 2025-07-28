clc
clear all
close all







%% 场群风资源输入设置
% 风向扇区
wind_sectors=8;

%风向划分，正北-东北：0-45度；东北-东：45-90度；东-东南：90-135度；东南-南：135-180；南-西南：180-225度；西南-西：225-270度；西-西北：270-325度；西北-北：325-360度
%1.正北-东北：0-45度
windfarmcluster_wind_direction_N_NE=0:5:45;
windfarmcluster_wind_speed_N_NE=3:1:12;

%2.东北-东：45-90度
windfarmcluster_wind_direction_NE_E=45:5:90;
windfarmcluster_wind_speed_NE_E=3:1:12;


%3.东-东南：90-135度
windfarmcluster_wind_direction_E_SE=90:5:135;
windfarmcluster_wind_speed_E_SE=3:1:10;

%4.东南-南：135-180
windfarmcluster_wind_direction_SE_S=135:5:180;
windfarmcluster_wind_speed_SE_S=3:1:10;

%5.南-西南：180-225度
windfarmcluster_wind_direction_S_SW=180:5:225;
windfarmcluster_wind_speed_S_SW=3:1:12;

%6.西南-西：225-270度
windfarmcluster_wind_direction_SW_W=225:5:270;
windfarmcluster_wind_speed_SW_W=3:1:5;


%7.西-西北：270-325度
windfarmcluster_wind_direction_W_NW=270:5:325;
windfarmcluster_wind_speed_W_NW=3:1:5;

%8.西北-北：325-360度
windfarmcluster_wind_direction_NW_N=325:5:360;
windfarmcluster_wind_speed_NW_N=3:1:5;




%% 初始化偏航设置
matrix=zeros(1,159);


%% 不同类型机组的设置
sqz_1=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','WindField','Range','B10:B10'));              %青州1风机个数
sqz_2=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','WindField','Range','B11:B11'));              %青州2风机个数
sqz_3=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','WindField','Range','B12:B12'));              %青州3风机个数
turbine_diameter_vector=[cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B3:B3')),...
                cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','E3:E3')),...
                cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','H3:H3'))];
turbine_hub_height_vector=[cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B4:B4')),...
                cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','E4:E4')),...
                cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','H4:H4'))];
rated_power_vector=[cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B13:B13')),...
                cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','E13:E13')),...
                cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','H13:H13'))];
life_total_vector=[cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B15:B15')),...
                cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','E15:E15')),...
                cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','H15:H15'))];
repair_c_vector=[cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B16:B16')),...
                cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','E16:E16')),...
                cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','H16:H16'))];
sqz_12=sqz_1+sqz_2;


without_optimization_result=struct();
optimized_resulut=struct();

% for i=1:length(windfarmcluster_wind_direction_N_NE)
%     for j=1:length(windfarmcluster_wind_speed_N_NE)
%% 初始化场群
swi=SmartWindInterface_yaw(sqz_12,turbine_diameter_vector,turbine_hub_height_vector,rated_power_vector,life_total_vector,repair_c_vector,matrix,0,3.5);
rng("default")
% swi.windfield.wind_direction=270;
swi.windfield.wake.velocity_model='Huadian';
swi.windfield.wake.deflection_model='Huadian';
swi.windfield.wake.turbulence_model='Huadian';
swi.windfield.enable_wfr='No';




% yaw_angles=30*rand(1,159);
% 
% 
% % 先热机
% swi.set_yaw_angles(yaw_angles);
% swi.calculate_wake();
% 
% % 第一次正式测
% t1 = tic;
% swi.set_yaw_angles(yaw_angles);
% swi.calculate_wake();
% t_elapsed1 = toc(t1);
% 
% % 第二次正式测
% t2 = tic;
% swi.set_yaw_angles(yaw_angles);
% swi.calculate_wake();
% swi.get_farm_power();
% t_elapsed2 = toc(t2);
% 
% fprintf('Run1 = %.4f s, Run2 = %.4f s (差值 %.2f%%)\\n', ...
%         t_elapsed1, t_elapsed2, 100*(t_elapsed1-t_elapsed2)/t_elapsed1);







% 优化并计时（整体）
% tOuter = tic;
tic
swi.yaw_optimization_gb();
toc

% % 设置 & 使用最大进程池
% pc = parcluster('local');
% delete(gcp('nocreate'));
% parpool(pc, pc.NumWorkers);   % 最大进程数
% % 或者：parpool('threads');


end_time=tic;
indexes=1:1:159;
a=swi.windfield.not_affecting_turbines;
b=setdiff(indexes, a); 
yaw_angles=rand(length(b),1);
[power, p_grad]=swi.cost_function_withgrad_block_parallel(swi.windfield.uwake_without_affect,yaw_angles,swi.windfield.wake_matrix_add, b);
end_time_total=toc(end_time);

fprintf('计算一次梯度: %.3f s\n', end_time_total);






% % elapsedOuter = toc(tOuter);

% 
% % 后续结果获取
% p2 = swi.get_farm_power();
% y2 = swi.get_yaw_angles();


















% % fprintf('外部脚本总调用耗时: %.3f s\n', elapsedOuter);











