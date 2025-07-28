classdef Turbine < handle
    
    properties
        rotor_diameter
        hub_height %=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B4:B4'));
        rated_power
        life_total
        repair_c
        
        %%%% 整个风场的环境湍流强度
        wind_field_turbulence_intensity=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Windfield','Range','B3:B3'));

        %%%优化相关量
        optimization_period=5;                                                                                  %5分钟优化一次
        annual_average_power=5*10^8;                                                                            %机组的年平均发电功率
        past_comprehensive_fatigue_coefficient=0.9800;                                                          %优化前的综合疲劳系数
        %%%%%优化相关量


        delta_t=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B14:B14'))                %优化周期                                                                                           %优化周期
        blade_count=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B5:B5'))

        
        pP=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B6:B6'))                       %青州123
        pT=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B7:B7'))                       %青州123

        generator_efficiency=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B8:B8'))
        power_thrust_table_11
        power_thrust_table_68
        power_thrust_table_83

        yaw_angle=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B10:B10'))              %偏航角度
        axial_induction_factor=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B12:B12'));%%轴向诱导因子

        tilt_angle=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B11:B11'))             %倾斜角度

        tsr=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range','B9:B9'))                      %叶尖速比(青州123)

        consumed_comprehensive_fatigue_coefficient=0.0292                                                       %已运行的寿命损耗，假设各风机相同

        air_density=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','WindField','Range','B8:B8'))            %空气密度
        velocities_u=[]
        turbulence_ambient=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Windfield','Range','B5:B5'));
        turbulence=[];
        efficiency_value=[]
        rotor_radius
        grid
    end
    
    properties (SetAccess=immutable,Hidden)
        ending_table_11=75;
        ending_table_68=135;
        ending_table_83=196;
    end  
    
    properties (Constant)
        points_turbine_grid=1;%风机(转子)网格点
    end

    properties (Dependent)       
        %grid
        %grid_with_yaw
        average_velocity
        Cp
        Ct
        aI%初始
        power
        comprehensive_fatigue_coefficient
        single_turbine_objective
        
    end

    methods
        %% 获取当前风机的转子直径和轮毂高度
        function obj=Turbine(sqz_12,turbine_diameter_vector,turbine_hub_height_vector,rated_power_vector,life_total_vector,repair_c_vector,count_tn)
            s_qz_12=sqz_12; 
            qz_3_count=table2array(readtable('青州3风机布局表.xlsx'));
            
            t_qz_12=1:1:s_qz_12;
            t_qz_3_1=92+find(qz_3_count(:,4)==6.8);
            t_qz_3_2=92+find(qz_3_count(:,4)==8.3);


        %% 确定风机类型，并赋值相应参数
            if ismember(count_tn,t_qz_12)
                obj.rotor_diameter=turbine_diameter_vector(1);
                obj.hub_height=turbine_hub_height_vector(1);
                obj.rated_power=rated_power_vector(1);
                obj.life_total=life_total_vector(1);
                obj.repair_c=repair_c_vector(1);
            elseif ismember(count_tn,t_qz_3_1)
                obj.rotor_diameter=turbine_diameter_vector(2);
                obj.hub_height=turbine_hub_height_vector(2);
                obj.rated_power=rated_power_vector(2);
                obj.life_total=life_total_vector(2);
                obj.repair_c=repair_c_vector(2);
            elseif ismember(count_tn,t_qz_3_2)
                obj.rotor_diameter=turbine_diameter_vector(3);
                obj.hub_height=turbine_hub_height_vector(3);
                obj.rated_power=rated_power_vector(3);
                obj.life_total=life_total_vector(3);
                obj.repair_c=repair_c_vector(3);
            end

       
        % % category_map = zeros(1, max_count_tn);
        % % category_map(t_qz_12) = 1;
        % % category_map(t_qz_3_1) = 2;
        % % category_map(t_qz_3_2) = 3;
        % % 
        % % category = category_map(count_tn);
        % % 
        % % switch category
        % %     case 1
        % %         idx = 1;
        % %     case 2
        % %         idx = 2;
        % %     case 3
        % %         idx = 3;
        % %     otherwise
        % %         error('Unrecognized count_tn');
        % % end
        % % 
        % % obj.rotor_diameter = turbine_diameter_vector(idx);
        % % obj.hub_height     = turbine_hub_height_vector(idx);
        % % obj.rated_power    = rated_power_vector(idx);
        % % obj.life_total     = life_total_vector(idx);
        % % obj.repair_c       = repair_c_vector(idx);



        %% Cp-CT表(不同型号机组对应不同表)
            obj.power_thrust_table_11=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range',sprintf('A18:C%d',obj.ending_table_11)));%功率系数表
            obj.power_thrust_table_68=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range',sprintf('A78:C%d',obj.ending_table_68)));%功率系数表
            obj.power_thrust_table_83=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Turbine','Range',sprintf('A139:C%d',obj.ending_table_83)));%功率系数表
        end

        %%%%%% %%%%%
        %% 重新定义的转子直径-更新转子半径
        function set.rotor_diameter(obj,d)
              obj.rotor_diameter=d;
              obj.update_radius;
        end

        %% 重新定义的转子半径+网格点的坐标，转子半径变化，风机转子处计速网格点坐标变化
        function set.rotor_radius(obj,r)
              obj.rotor_radius=r;
              obj.update_grid;
        end


%%%%% 功率推力系数查表 %%%%% 
%% 确定当前风速下对应的功率系数，基于表数据的线性插值
        function cp = fCp(obj,at_wind_speed)
            switch obj.rated_power
                case 11
                    power_thrust_table = obj.power_thrust_table_11;
                case 6.8
                    power_thrust_table = obj.power_thrust_table_68;
                case 8.3
                    power_thrust_table = obj.power_thrust_table_83;
            end

            cp_column=power_thrust_table(:,2);%提取不同风速下对应的功率系数
            wind_speed_column=power_thrust_table(:,1);%提取风速
            if at_wind_speed < wind_speed_column(1) || at_wind_speed > wind_speed_column(end)
                cp=0;
            else 
                cp=interp1(wind_speed_column,cp_column,at_wind_speed,'linear');
            end
        end

%% 确定当前风速下对应的推力系数，基于表数据的线性插值，功率为0，推力系数不为0
        function ct = fCt(obj,at_wind_speed)
            switch obj.rated_power
                case 11
                    power_thrust_table = obj.power_thrust_table_11;
                case 6.8
                    power_thrust_table = obj.power_thrust_table_68;
                case 8.3
                    power_thrust_table = obj.power_thrust_table_83;
            end
            ct_column=power_thrust_table(:,3);
            wind_speed_column=power_thrust_table(:,1);
            if at_wind_speed < wind_speed_column(1) 
                ct=0.99;
            elseif at_wind_speed > wind_speed_column(end)
                ct=0.0001;
            else
                ct=interp1(wind_speed_column,ct_column,at_wind_speed,'linear');
            end
        end
%%%%功率推力系数查表%%%%%%%%%





%%%%%%%%更新当前风机关联点处的风速和附加湍流强度%%%%%%%%%
          %% 确定当前风机所有相关网格点处的风速
          function data=calculate_turbine_velocities(obj,local_wind_speed,coord,x,y,z,windfield,nt)
            y_rel=obj.grid(:,1);                                        %所有网格点的y坐标，从左到右
            z_rel=obj.grid(:,2);                                        %所有网格点的x坐标，从上到下          
            if strcmp(windfield.enable_wfr,'No')
                data=zeros(length(y_rel),1);                            %size:16x1
                selected_layer=local_wind_speed(:,:,nt);                %nt表示的是第nt台风机，local_wind_speed表示各台风机划分的网格点处的风速
                for i=1:length(y_rel)
                   data(i)=selected_layer(obj.grid(i,3));               %local_wind_speed
                end
            else
            distance=zeros(numel(x),length(y_rel));                     %numel(x)代表了全部风机的网格点个数如16x4，size:(16x4)x16
            idx=zeros(length(y_rel),1);                                 %当前风机所有网格点的个数size:16x1
            data=zeros(length(y_rel),1);                                %size:16x1
                for i=1:length(y_rel)                                   
                    for j=1:numel(x)                                    %全部风机的所有网格点
                        distance(j,i)=sqrt((coord(1)-x(j))^2+(coord(2)+y_rel(i)-y(j))^2+(coord(3)+z_rel(i)-z(j))^2);
                    end
                    [~,idx(i)]=min(distance(:,i));
                    data(i)=local_wind_speed(idx(i));
                end
            end
          end

          %% 确定当前风机所有相关网格点处的附加湍流强度
          function data=calculate_turbine_addedturbulenceintensity(obj,local_turbulence,coord,x,y,z,windfield,nt)
            y_rel=obj.grid(:,1);                                                                           %所有网格点的y坐标，从左到右
            z_rel=obj.grid(:,2);                                                                           %所有网格点的x坐标，从上到下          
            if strcmp(windfield.enable_wfr,'No')
                data=zeros(length(y_rel),1);%size:16x1
                selected_layer=local_turbulence(:,:,nt);%nt表示的是第nt台风机，local_wind_speed表示各台风机划分的网格点处的风速
                for i=1:length(y_rel)
                   data(i)=selected_layer(obj.grid(i,3));%local_wind_speed
                end
            else
            distance=zeros(numel(x),length(y_rel));%numel(x)代表了全部风机的网格点个数如16x4，size:(16x4)x16
            idx=zeros(length(y_rel),1);%当前风机所有网格点的个数size:16x1
            data=zeros(length(y_rel),1);%size:16x1
                for i=1:length(y_rel)%
                    for j=1:numel(x)%全部风机的所有网格点
                        distance(j,i)=sqrt((coord(1)-x(j))^2+(coord(2)+y_rel(i)-y(j))^2+(coord(3)+z_rel(i)-z(j))^2);
                    end
                    [~,idx(i)]=min(distance(:,i));
                    data(i)=local_turbulence(idx(i));
                end
            end
          end
%%%%%%%%%更新当前风机关联点处的风速和附加湍流强度%%%%%%%%%




%%%%%%%%更新当前风机关联点处的风速和湍流强度%%%%%%%%%%%%
        %% 更新当前风机所有网格点处的风速(考虑尾流影响后)
        function obj=update_velocities(obj,u_wake,windfield,coord,rotated_x,rotated_y,rotated_z,u_initial,nt)
            local_wind_speed_u=u_initial-u_wake;
            obj.velocities_u=calculate_turbine_velocities(obj,local_wind_speed_u,coord,rotated_x,rotated_y,rotated_z,windfield,nt);
        end

        % %% 更新当前风机所有网格点处的湍流强度(考虑尾流影响后)
        % function obj=update_turbulences(obj,u_wake,windfield,coord,rotated_x,rotated_y,rotated_z,u_initial,nt)
        %     local_wind_speed_u=u_initial-u_wake;
        %     obj.velocities_u=calculate_turbine_velocities(obj,local_wind_speed_u,coord,rotated_x,rotated_y,rotated_z,windfield,nt);
        % end
%%%%%%%%更新当前风机关联点处的风速和湍流强度%%%%%%%%%%%%

        
        %%  更新当前风机所有网格点处的平均风速
        function average_velocity=get.average_velocity(obj)
            average_velocity=(mean((obj.velocities_u).^3))^(1/3);%%开3次根号
        end
        %% 更新功率系数
        function Cp=get.Cp(obj)
            % pW=obj.pP/3;
            yaw_effective_velocity=obj.average_velocity*(cosd(obj.yaw_angle));
            Cp=fCp(obj,yaw_effective_velocity);

            
            % aif=obj.axial_induction_factor;
            % Cp=4*aif*(1-aif)^2*(cosd(obj.yaw_angle))^obj.pP;


        end
        %% 更新推力系数
        function Ct=get.Ct(obj)
            yaw_effective_velocity=obj.average_velocity*(cosd(obj.yaw_angle));
            Ct=fCt(obj,yaw_effective_velocity);
            
           % aif=obj.axial_induction_factor;
           %  Ct=4*aif*(1-aif)*(cosd(obj.yaw_angle))^obj.pT;


        end
        %% 湍流计算过程量
        function aI=get.aI(obj)
            aI=0.5/cosd(obj.yaw_angle)*(1-sqrt(1-obj.Ct*cosd(obj.yaw_angle)));
        end

        %% 获取风机的功率        
        function power=get.power(obj)
            %计算效能
            % eff_value=obj.efficiency_value;
            % eff_value=0.98;
            
            %pW=obj.pP/3;
            yaw_effective_velocity=obj.average_velocity*(cosd(obj.yaw_angle));
           
            cptmp=obj.Cp;
            
            power=0.5*obj.air_density*pi*obj.rotor_radius^2*obj.generator_efficiency*yaw_effective_velocity^3*cptmp;
        end

        %% 获取风机累积的疲劳损伤
        function comprehensive_fatigue_coefficient=get.comprehensive_fatigue_coefficient(obj)
            % 累计疲劳系数
            accumulated_c=obj.consumed_comprehensive_fatigue_coefficient;

            % 待优化的疲劳系数
            optimized_c_power=obj.power*obj.delta_t/(obj.rated_power*obj.life_total*(1+obj.repair_c));
            optimized_c_turbulence=1.5*obj.turbulence*obj.delta_t/(obj.rated_power*obj.life_total*(1+obj.repair_c));
            optimized_c=optimized_c_power+optimized_c_turbulence;

            % 剩余寿命比
            comprehensive_fatigue_coefficient=1-(accumulated_c+optimized_c);
        end  

        %% 单机优化目标
        function single_turbine_objective=get.single_turbine_objective(obj)
            %%尽可能最大化1)正常所有风机的发电量+2)由各机组寿命优化折算得到的发电量的总和
                 single_turbine_objective=obj.power*obj.optimization_period+obj.annual_average_power*(obj.comprehensive_fatigue_coefficient-obj.past_comprehensive_fatigue_coefficient)*obj.optimization_period;

        end

        %% 获取风机的效能
        function obj=update_efficiency_value(obj,k_to_m_vector,m_to_e_vector)
         %动能到机械能&机械能到电能转化的系数
          kinetic_to_mechanichal_coefficient=[0.0424, 0.042, 0.0425, 0.0423, 0.0424, 0.0417, 0.0425, 0.0426, 0.0425, 0.0412, 0.0426, 0.0426, 0.042, 0.0308, 0.0424, 0.0426, 0.0427, 0.0426, 0.042, 0.0421, 0.0387, 0.0424, 0.0424, 0.0422];
          mechanical_to_electrical_coefficient=[0.142,0.143,0.1429,0.1427,0.1427,0.1426,0.144];

         %动能到机械能&机械能到电能效能计算
         efficiency_k_to_m=kinetic_to_mechanichal_coefficient'.*k_to_m_vector;
         efficiency_m_to_e=mechanical_to_electrical_coefficient'.*m_to_e_vector;

         %总的效能值
         obj.efficiency_value=efficiency_k_to_m*efficiency_m_to_e;
        
        end


        %% 更新风机的湍流强度
        % function obj=update_turbulence_intensity(obj,windfield,wake,turbine_coord,wake_coord,turbine_wake,area)
        %     i=wake.ti_initial;
        %     constant=wake.ti_constant;
        %     ai=wake.ti_ai;
        %     downstream=wake.ti_downstream;
        %     added_ti=area*constant*(turbine_wake.aI^ai)*(windfield.turbulence_intensity^i)*((turbine_coord(1)-wake_coord(1))/obj.rotor_diameter)^(downstream);
        %     obj.turbulence=sqrt(obj.turbulence^2+added_ti^2);
        % end

        %% 更新当前风机i的入流湍流强度
        function obj=update_turbulence_intensity(obj,windfield,coord,rotated_x,rotated_y,rotated_z,u_turbulence_wake,nt)
            local_turbulence=u_turbulence_wake;
            
            data=calculate_turbine_velocities(obj,local_turbulence,coord,rotated_x,rotated_y,rotated_z,windfield,nt);
            
            added_ti=sqrt(sum(data.^2));
            % added_ti=area*constant*(turbine_wake.aI^ai)*(windfield.turbulence_intensity^i)*((turbine_coord(1)-wake_coord(1))/obj.rotor_diameter)^(downstream);
            obj.turbulence=sqrt(obj.turbulence_ambient^2+added_ti^2);
        end
        
        function obj=update_radius(obj)
            obj.rotor_radius=obj.rotor_diameter/2;
        end

        %% 旋转后更新各风机处所有网格点的坐标
        % function obj=update_grid(obj)                                                               %从上到下，从左到右把所有网格点的坐标列出来
        %     num_points=sqrt(obj.points_turbine_grid);                                               %网格水平、垂直方向划分点的个数
        %     horizontal=linspace(-obj.rotor_radius,obj.rotor_radius,num_points);                     %网格水平位置坐标划分
        %     vertical=linspace(obj.rotor_radius,-obj.rotor_radius,num_points);                       %网格垂直位置坐标划分
        %     obj.grid=zeros(obj.points_turbine_grid,3);
        %     count_points=linspace(1,obj.points_turbine_grid,obj.points_turbine_grid);
        %     obj.grid(:,3)=count_points;                                                             %grid的第三列为网格点的索引
        %     for v=1:length(vertical)
        %         for h=1:length(horizontal)
        %              obj.grid(num_points*(v-1)+h,1)=[horizontal(h)];
        %              obj.grid(num_points*(v-1)+h,2)=[vertical(v)];
        %         end
        %     end
        %     for n=obj.points_turbine_grid:-1:1
        %         if sqrt(obj.grid(n,1)^2+obj.grid(n,2)^2)>obj.rotor_radius
        %         obj.grid(n,:)=[];
        %         end
        %     end
        % end 

        %% 风机转子平面处的网格点：取法一只在转子中心横轴上取 3 个点（左、中、右）
        function obj = update_grid(obj)

            num_points   = obj.points_turbine_grid;                             
            xs           = linspace(-obj.rotor_radius, ...
                                     obj.rotor_radius, ...
                                     obj.points_turbine_grid);         % 三个水平位置[-R, 0, R]

            obj.grid     = zeros(num_points, 3);          % [x, y, index]
            obj.grid(:,1)= xs;                            % 横向坐标
            obj.grid(:,2)= 0;                             % 纵向（垂直）坐标，全为 0（中心线）
            obj.grid(:,3)= (1:num_points)';               
        end 
    end
end


