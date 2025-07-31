clc
clear
close all

rng("default");
% 4.07361843196590	28.5289596853781	48.6349340814675	76.5668792806951	
% 99.1617962311270	120.487702024997	145.392491094335	170.734407596025	
% 196.787534177172	220.824442675996	240.788065408388	268.852963908803	
% 292.785834741215	314.426878243614	340.001402344444
%% NOTE:
% wd(1)=15.2 + 0.1*rand pass check
% wd(6)=115 + 0.1*rand pass check
% wd(8)=158.8 + 0.1*rand pass check
% wd(9)=203.7 + 0.1*rand pass check
% wd(13)=312.7 + 0.1*rand pass check
% wd(15)=345.1 + 0.1*rand pass check


% 3	 5.33333333333333	7.66666666666667	10	12.3333333333333
% 14.6666666666667	17	19.3333333333333	21.6666666666667	24
%% 定义测试集
wind_sectors = 15; % from 0 to 360 degrees
velocity_sectors = 10; % from 3m/s to 25m/s
wind_direction = linspace(0, 360*(1-1/wind_sectors), wind_sectors);
wind_direction = wind_direction + 5*rand(size(wind_direction));
vel_sets = linspace(3, 24, velocity_sectors);
qingzhou12_agc_perc = [0.5 1]; % 50%与100%基准功率
qingzhou3_agc_perc = [0.5 1];

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

%% 测试数据记录
objectives = zeros(length(vel_sets), 4); % wind_direction, vel_sets, power_agc(4 types)
times = zeros(length(vel_sets), 4); % 记录每次计算的时间
err_list = zeros(1, length(vel_sets)); % 错误列表
%% 运行测试
% 检查是否有并行池，如有则直接使用，否则创建一个新的并行池
if isempty(gcp('nocreate'))
    parpool('local', 24); % 启动并行池
else
    fprintf('Using existing parallel pool.\n');
end
restore_mode = 1;
wd_chk_point = 345.1; % 检查点，避免重复计算
for wd = 1
    for vel = 1:length(vel_sets)
        fprintf('Wind direction: %.1f, Wind speed: %.1f\n', wd_chk_point(wd), vel_sets(vel));
        swi = SmartWindInterface_yaw(sqz_12, turbine_diameter_vector, turbine_hub_height_vector, rated_power_vector,...
         life_total_vector, repair_c_vector, matrix, wd_chk_point(wd), vel_sets(vel));
        swi.windfield.wake.velocity_model = 'Huadian';
        swi.windfield.wake.deflection_model = 'Huadian';
        swi.windfield.wake.turbulence_model = 'Huadian';
        swi.windfield.enable_wfr = 'No';
        try
            swi.calculate_wake();
            power_baseline = swi.get_farm_power();
            qingzhou12_baseline = swi.get_farm_qingzhou12_power();
            qingzhou3_baseline = swi.get_farm_qingzhou3_power();
            qingzhou12_agc = qingzhou12_baseline .* qingzhou12_agc_perc; % 50%与100%基准功率
            qingzhou3_agc = qingzhou3_baseline .* qingzhou3_agc_perc;
            tic
            swi.yaw_optimization_tracking_life_ipopt(qingzhou12_agc(1), qingzhou3_agc(1));
            times((wd - 1) * length(vel_sets) + vel, 1) = toc;
            objectives((wd - 1) * length(vel_sets) + vel, 1) = swi.get_farm_power();
            disp('pass test #1')
            % tic
            % swi.yaw_optimization_tracking_life_ipopt(qingzhou12_agc(1), qingzhou3_agc(2));
            % times((wd - 1) * length(vel_sets) + vel, 2) = toc;
            % objectives((wd - 1) * length(vel_sets) + vel, 2) = swi.get_farm_power();
            % disp('pass test #2')
            % tic
            % swi.yaw_optimization_tracking_life_ipopt(qingzhou12_agc(2), qingzhou3_agc(1));
            % times((wd - 1) * length(vel_sets) + vel, 3) = toc;
            % objectives((wd - 1) * length(vel_sets) + vel, 3) = swi.get_farm_power();
            % disp('pass test #3')
            % tic
            % swi.yaw_optimization_tracking_life_ipopt(qingzhou12_agc(2), qingzhou3_agc(2));
            % times((wd - 1) * length(vel_sets) + vel, 4) = toc;
            % objectives((wd - 1) * length(vel_sets) + vel, 4) = swi.get_farm_power();
            % disp('pass test #4')
            % 测试通过
            fprintf('Wind direction: %.1f, Wind speed: %.1f - Test passed.\n', wd_chk_point(wd), vel_sets(vel));
        catch
            err_list((wd - 1) * length(vel_sets) + vel) = -1; % 记录错误
            % 测试不通过
            fprintf('Wind direction: %.1f, Wind speed: %.1f - Test failed.\n', wd_chk_point(wd), vel_sets(vel));
        end
    end
end
