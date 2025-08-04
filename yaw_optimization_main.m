clear
rng('default')
% 导入风场相关数据
load('Input/input_data.mat');

% 选择数据点
wd = 1; vel = 1;
% 初始化风场相关class
swi = SmartWindInterface_yaw(sqz_12, turbine_diameter_vector, turbine_hub_height_vector, rated_power_vector,...
    life_total_vector, repair_c_vector, matrix, wind_direction(wd), vel_sets(vel));
swi.windfield.wake.velocity_model = 'Huadian';
swi.windfield.wake.deflection_model = 'Huadian';
swi.windfield.wake.turbulence_model = 'Huadian';
swi.windfield.enable_wfr = 'No';
swi.calculate_wake();

% 计算目标功率
power_baseline = swi.get_farm_power();
qingzhou12_baseline = swi.get_farm_qingzhou12_power();
qingzhou3_baseline = swi.get_farm_qingzhou3_power();
qingzhou12_agc = qingzhou12_baseline .* qingzhou12_agc_perc; % 50%与100%基准功率
qingzhou3_agc = qingzhou3_baseline .* qingzhou3_agc_perc;

% 检查是否有并行池，如有则直接使用，否则创建一个新的并行池
if isempty(gcp('nocreate'))
    parpool('local', 32); % 启动并行池
else
    fprintf('Using existing parallel pool.\n');
end

% 打印输入：风向，风速，青州12与青州3的AGC功率百分比
fprintf('本次执行结果为：风向: %.0f, 风速: %.1f, 青州12风场AGC: %.2f%%, 青州3风场AGC: %.2f%%\n', ...
    wind_direction(wd), vel_sets(vel), 100*qingzhou12_agc_perc(1), 100*qingzhou3_agc_perc(1));
    
%% 计时开始
tic
swi.yaw_optimization_tracking_life_ipopt(qingzhou12_agc(1), qingzhou3_agc(1));
time = toc; % 计时结束并记录
objectives((wd - 1) * length(vel_sets) + vel, 1) = swi.get_farm_power();
% 打印运行时间
fprintf('本次执行周期为: %.2fs\n', time);
% 判定是否小于29s
if times((wd - 1) * length(vel_sets) + vel, 1) <= 29
    fprintf('本次测试周期时间小于29秒\n');
end