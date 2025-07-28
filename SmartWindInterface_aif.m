classdef SmartWindInterface_aif <handle

    properties
        layout_x
        layout_y
        turbine_status%风机状态
        windfield%风场对象
        imaging
        energy
        aif_lower=0;%%最小的偏航动作角
        aif_upper=1/3;%%最大的偏航动作角
        ya_options=10
        
        


    end

    properties (SetAccess=immutable,Hidden)
        ending_row=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Layout',...  %定义ending_row
            'Range','G5:G5'));
    end


    methods

        function obj=SmartWindInterface_aif(turbine_diameter_vector,turbine_hub_height_vector,rated_power_vector,life_total_vector,repair_c_vector,yaw_matrix)
            obj.layout_x=cell2mat(readcell('inputs_all_fields.xlsx','Sheet',...
                'Layout','Range',sprintf('B3:B%d',obj.ending_row)));               %场群中风机布局的x坐标
            obj.layout_y=cell2mat(readcell('inputs_all_fields.xlsx','Sheet'...
                ,'Layout','Range',sprintf('C3:C%d',obj.ending_row)));              %场群中风机布局的y坐标
            wake=Wake;
            n_turbines=length(obj.layout_x);                                       %场群中风机的个数
            turbines=cell(n_turbines,1);
            for i=1:n_turbines
                turbines{i,1}=Turbine(turbine_diameter_vector,turbine_hub_height_vector,rated_power_vector,life_total_vector,repair_c_vector,i); 
                turbines{i,1}.yaw_angle=yaw_matrix(i);                             %风机对象类
            end
            turbinechart=Turbinechart(obj.layout_x,obj.layout_y,turbines);         %Turbinechart对象类
            obj.turbine_status=ones(n_turbines,1);                                 %风机状态
            obj.windfield=Windfield(turbinechart,wake);                            %场群对象类
            obj.imaging=Imaging;
            obj.energy=Energy;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%SIMULATION_PART%%%%%%%%%%%%%%%%%%%%%%%%%%%

        




        %SET_LAYOUT is a function to set a layout different from the one 
        %initialized in the excel file. The input must be a matrix nx2 with
        %n the number of turbines; for each turbine coordinate x and y must
        %be expressed. Alternatively, use the layout generators to set an 
        %array of turbines (eventually rotated) or a random layout.
        %Remember to set the layout before setting any other attribute of
        %other classes (velocity,direction,...), otherwise any other
        %property will be reset.

        %% 自定义风场布局
        function obj=set_layout(obj,matrix)
            obj.layout_x=matrix(:,1);
            obj.layout_y=matrix(:,2);
            wake=Wake;
            n_turbines=length(obj.layout_x);
            turbines=cell(n_turbines,1);
            for i=1:n_turbines
                turbines{i,1}=Turbine;
            end
            turbinechart=Turbinechart(obj.layout_x,obj.layout_y,turbines);
            obj.turbine_status=ones(n_turbines,1);
            obj.windfield=Windfield(turbinechart,wake);
        end
        

        %EXCLUDE_TURBINES is a function to exclude some turbines temporary
        %from a layout and make them "transparent". This could be useful in
        %reality due to faulted turbines that cannot be involved in the 
        %calculation process and in the yaw optimization. The input must be
        %a vector of integers that correspond to the numbers of the
        %excluded turbines (in the same order as presented in layout).
        %Remember to call this function before setting any other attribute
        %of other classes. This function has the priority with respect to
        %other commands except for the function "set_layout".

        % 排除故障涡轮机(优先设置)
        function obj=exclude_turbines(obj,exc)
            obj.reset_farm_keep_layout();
            obj.turbine_status=ones(length(obj.layout_x),1);
            for i=1:length(exc)
                turb=exc(i);
                sz=length(obj.windfield.turbinechart.turbines...
                    {turb,1}.power_thrust_table);
                obj.windfield.turbinechart.turbines...
                    {turb,1}.power_thrust_table(:,2)=repelem(0.0001,sz);
                obj.windfield.turbinechart.turbines...
                    {turb,1}.power_thrust_table(:,3)=zeros(sz,1);
                obj.turbine_status(exc(i))=0;
            end
        end
        
        
        %CALCULATE_WAKE is the core function of the program. It enables to 
        %calculate the wind flow speed at the turbine points thanks to 
        %several wake models listed in the Excel file. As a result, also
        %the characteristic wind speeds of each turbine are calculated and
        %stored in the windfield object, together with the power of the
        %turbines, the turbulence of the turbines and the farm power. To
        %speed up the calculations (for example in case of plant online
        %optimization) it is recommended NOT to enable the resolution so
        %that the wind flow field is calculated only at the turbine
        %points. Instead, to have a graphical representation of the whole
        %field the resolution option must be enabled. Clearly this will
        %affect significantly the computational time, if the resolution is
        %low (e.g [50 30 10]) also the accuracy of the calculated 
        % parameters. On the other hand, [250 150 50] will be a medium 
        %resolution, while [500 300 100] will be a high resolution   




        
        %% 计算尾流
        function obj=calculate_wake(obj)
            obj=obj.windfield.calculatewake();
        end
        
        %CALCULATE_NOWAKE is a function mainly used to compute the flow
        % conditions in the irrealistic condition of absence of wake. This
        % is useful to compute the wake losses and the efficiency of the
        % wind farm. Clearly, at a given velocity the power produced by 
        % each turbine will be the same regardless the direction
        
        %% 无尾流非真实条件下的流动状况
        function obj=calculate_nowake(obj)
            obj=obj.windfield.calculatenowake();
        end
        

        %RESET_FARM is a function to restart the interface object so that
        %every modification done in the Matlab working file to the 
        % properties of the object are cancelled and the inputs are again
        % the same of the input file. 

        %% 接口对象重启函数，原来的对象属性修改清0
        function obj=reset_farm(obj)
            obj.layout_x=cell2mat(readcell('inputs_all_fields.xlsx','Sheet',...
                'Layout','Range',sprintf('B3:B%d',obj.ending_row)));
            obj.layout_y=cell2mat(readcell('inputs_all_fields.xlsx','Sheet',...
                'Layout','Range',sprintf('C3:C%d',obj.ending_row)));
            wake=Wake;
            n_turbines=length(obj.layout_x);
            turbines=cell(n_turbines,1);
            for i=1:n_turbines
                turbines{i,1}=Turbine;
            end
            turbinechart=Turbinechart(obj.layout_x,obj.layout_y,turbines);
            obj.turbine_status=ones(n_turbines,1);
            obj.windfield=Windfield(turbinechart,wake);
        end
        

        %RESET_FARM_KEEP_LAYOUT is a function to restart the interface
        %properties with the exception of the layout modified in the Matlab
        %working file. As a result, with the exception of the layout, all
        %other input properties will be the same of the input Excel file.

        %% 重启接口对象属性(除了Matlab工作文件中修改的布局除外)，除了布局以外，所有其他输入属性都重置为Excel文件记录
        function obj=reset_farm_keep_layout(obj)
            wake=Wake;
            n_turbines=length(obj.layout_x);
            turbines=cell(n_turbines,1);
            for i=1:n_turbines
                turbines{i,1}=Turbine;
            end
            turbinechart=Turbinechart(obj.layout_x,obj.layout_y,turbines);
            obj.turbine_status=ones(n_turbines,1);
            obj.windfield=Windfield(turbinechart,wake);
        end
        

        %GET_YAW_ANGLES is a useful function that outputs a vector with the
        %yaw angles of all the turbines, with the same order the turbines
        %are listed in the layout properties. This is a shortcut, since the
        %yaw angles are stored in different turbine objects and recalling
        %them would require a for loop in the code each time.

        % %% 存储所有涡轮机的偏航角
        % function yaw_vector=get_yaw_angles(obj)
        %     yaw_vector=zeros(1,length(obj.layout_x));
        %     for i=1:length(yaw_vector)
        %         yaw_vector(i)=obj.windfield.turbinechart.turbines...
        %             {i,1}.yaw_angle;
        %     end
        % end
        
        %% 存储所有涡轮机的轴向诱导因子
        function aif_vector=get_aif_angles(obj)
            aif_vector=zeros(1,length(obj.layout_x));
            for i=1:length(aif_vector)
                aif_vector(i)=obj.windfield.turbinechart.turbines...
                    {i,1}.axial_induction_factor;
            end
        end




        %GET_TILT_ANGLES is a function that has the same goal of the
        %'get_yaw_angles()' function. However, no tilt angle steering is
        %employed at the moment in this program, but this could be an
        %interesting option for the future of wind turbines.

        %% 存储所有风机的倾斜角度
        function tilt_vector=get_tilt_angles(obj)
            tilt_vector=zeros(1,length(obj.layout_x));
            for i=1:length(tilt_vector)
                tilt_vector(i)=obj.windfield.turbinechart.turbines...
                    {i,1}.tilt_angle;
            end
        end

        %SET_YAW_ANGLES is a useful function that, given an input vector of
        %yaw angles (one for each turbine), modifies the property of the
        %yaw angle for each turbine object.This is a necessary shortcut,
        %since otherwise a for loop would have been necessary to store the
        %values in each object. The yaw angles refer to each turbine with 
        %the same order turbines are listed in the layout properties

        % %% 存储所有涡轮机偏航角的输入向量(风机指定偏航角)
        % function obj=set_yaw_angles(obj,input_vector)
        %     if length(input_vector)~=length(obj.layout_x)
        %         error(['Turbines number and yaw angles number do' ...
        %             ' not correspond'])
        %     else
        %         for i=1:length(input_vector)
        %             obj.windfield.turbinechart.turbines...
        %                 {i,1}.yaw_angle=input_vector(i);
        %         end
        %     end
        % end

        %% 存储所有涡轮机轴向诱导因子的输入向量(风机指定轴向诱导因子)
        function obj=set_aif(obj,input_vector)
            if length(input_vector)~=length(obj.layout_x)
                error(['Turbines number and axial induction factor number do' ...
                    ' not correspond'])
            else
                for i=1:length(input_vector)
                    obj.windfield.turbinechart.turbines...
                        {i,1}.axial_induction_factor=input_vector(i);
                end
            end
        end




        %SET_TILT_ANGLES is a function that has the same goal of the
        %'set_yaw_angles()' function. However, no tilt angle steering is
        %employed at the moment in this program, but this could be an
        %interesting option for the future of wind turbines.

        %% 各风机指定的倾斜角
        function obj=set_tilt_angles(obj,input_vector)
            if length(input_vector)~=length(obj.layout_x)
                error(['Turbines number and tilt angles number do' ...
                    ' not correspond'])
            else
                for i=1:length(input_vector)
                    obj.windfield.turbinechart.turbines...
                        {i,1}.tilt_angle=input_vector(i);
                end
            end
        end
        
        
        %GET_TURBINES_POWER is a function that outputs a vector with the
        %power of each turbine. This function has to be called after
        %the wake calculation or yaw optimization to output the turbine
        %powers in a specific situation.The powers refer to each 
        %turbine with the same order turbines are listed in the layout
        %properties

        %% 存储各风机的功率
        function power_cell=get_turbines_power(obj)
            power_cell=zeros(length(obj.layout_x),1);
            for i=1:length(obj.layout_x)
                power_cell(i)=obj.windfield.turbinechart.turbines...
                    {i,1}.power;
            end
        end
        
        %GET_FARM_POWER is a function that outputs the total power produced
        %by the wind farm, obtained by summing all the individual powers
        %produced by each turbine. It has to be called after the wake
        %calculation or yaw optimization.

        %% 存储整个风场的功率
        function total_power=get_farm_power(obj)
            power_cell=zeros(length(obj.layout_x),1);
            for i=1:length(obj.layout_x)
                power_cell(i)=obj.windfield.turbinechart.turbines...
                    {i,1}.power;
            end
            total_power=sum(power_cell);
        end


         %% 存储各风机的寿命
        function life_cell=get_turbines_life(obj)
            life_cell=zeros(length(obj.layout_x),1);
            for i=1:length(obj.layout_x)
                life_cell(i)=obj.windfield.turbinechart.turbines...
                    {i,1}.left_life;
            end
        end
        
        %GET_FARM_POWER is a function that outputs the total power produced
        %by the wind farm, obtained by summing all the individual powers
        %produced by each turbine. It has to be called after the wake
        %calculation or yaw optimization.

        %% 存储整个风场的平均寿命
        function total_life=get_farm_life(obj)
            life_cell=zeros(length(obj.layout_x),1);
            for i=1:length(obj.layout_x)
                life_cell(i)=obj.windfield.turbinechart.turbines...
                    {i,1}.power;
            end
            total_life=sum(life_cell)/length(obj.layout_x);
        end


        %% 存储各风机的寿命+功率
        function objective_cell=get_turbines_objective(obj)
            objective_cell=zeros(length(obj.layout_x),1);
            for i=1:length(obj.layout_x)
                objective_cell(i)=obj.windfield.turbinechart.turbines...
                    {i,1}.single_turbine_objective;
            end
        end
        
        %GET_FARM_POWER is a function that outputs the total power produced
        %by the wind farm, obtained by summing all the individual powers
        %produced by each turbine. It has to be called after the wake
        %calculation or yaw optimization.

        %% 存储整个风场的功率+寿命
        function total_objective=get_farm_objective(obj)
            objective_power_cell=zeros(length(obj.layout_x),1);
            objective_life_cell=zeros(length(obj.layout_x),1);
            for i=1:length(obj.layout_x)
                objective_power_cell(i)=obj.windfield.turbinechart.turbines...
                    {i,1}.power;
                objective_life_cell(i)=obj.windfield.turbinechart.turbines...
                    {i,1}.left_life;

            end
            total_power_objective=sum(objective_power_cell);
            total_life_objective=mean(objective_life_cell)-sqrt(mean((objective_life_cell-mean(objective_life_cell)).^2));
            total_objective=total_power_objective+total_life_objective;
        end
        








        
        %GET_TURBINES_VELOCITY is a function that outputs a vector
        %containing as element the velocities at each turbine. Since in
        %this software the turbine are represented by several grid points
        %placed on the rotor surface, the cubic mean is calculated for each
        %turbine to get a representative speed for each turbine. It has 
        % to be called after the wake calculation or yaw optimization. 
        % The velocities refer to each turbine with the same order turbines
        % are listed in the layout properties

        %% 每台涡轮机的代表风速(由于网格的存在)
        function velocities_cell=get_turbines_velocity(obj)
            velocities_cell=zeros(length(obj.layout_x),1);
            for i=1:length(obj.layout_x)
                velocities_cell(i)=obj.windfield.turbinechart.turbines...
                    {i,1}.average_velocity;
            end
        end

        %% 每台涡轮机的湍流
        %GET_TURBINES_TURBULENCE is a function that outputs a vector
        %containing the turbulence value of each turbine calculated with 
        %Crespo-Hernandez model. Since some models (Jensen, Jimenez,
        %Multizone) do not require turbulence calculation, the command to
        %calculate turbulence has to be enabled manually. On the other hand
        %Gaussian models require the turbulence for calculations so the
        %turbulence values are automatically calculated. Clearly, this 
        %function has to be called after wake calculation or yaw
        %optimization. The velocities refer to each turbine with the same 
        %order turbinesare listed in the layout properties
        function turbulence_cell=get_turbines_turbulence(obj)
            turbulence_cell=zeros(length(obj.layout_x),1);
            for i=1:length(obj.layout_x)
                turbulence_cell(i)=obj.windfield.turbinechart.turbines...
                    {i,1}.turbulence;
            end
        end
        
        %SHOW_HORPLANE is a function that outputs a rendering of the wind
        %farm from above, cutting a plane at some z-axis value, that has to
        %be specified in the input. The most representative cutpoint is
        %obviously the hub height of the turbines. If the cutpoint is the
        %hub height and all the turbines have the same height also the
        %turbine sections are shown in the picture. If the turbines are red
        %, it means that they have been excluded. This function has to be
        %called after 'calculate_wake()' and resolution must be enabled.

        %% 轮毂高度处的风场切片，俯视图，涡轮机红色表示其被排除在外
        function pseudocolor=show_horplane(obj,cutpoint)
             hub_height_cell=zeros(length(obj.layout_x),1);
             for i=1:length(obj.layout_x)
                hub_height_cell(i)=obj.windfield.turbinechart.turbines...
                    {i,1}.hub_height;
             end
             pseudocolor=obj.imaging.z_view(obj.windfield,cutpoint);
             if all(hub_height_cell==cutpoint)
                hold on
                obj.imaging.plot_turbines(obj.layout_x,obj.layout_y,...
                    obj.windfield,obj);
             end
        end
        


        %SHOW_VERPLANE is a function that renders the points of the flow
        %field in the y-z plane, cutting it at some x-axis value, that has
        %to be specified in the input. This function has to be
        %called after 'calculate_wake()' and resolution must be enabled.
        
        %% y-z平面切片
        function pseudocolor=show_verplane(obj,cutpoint)
             pseudocolor=obj.imaging.x_view(obj.windfield,cutpoint);
        end

        %SHOW_CROSSPLANE is a function that renders the points of the flow
        %field in the x-z plane, cutting it at some y-axis value, that has
        %to be specified in the input. This function has to be
        %called after 'calculate_wake()' and resolution must be enabled.

        %% x-z平面切片
        function pseudocolor=show_crossplane(obj,cutpoint)
             pseudocolor=obj.imaging.y_view(obj.windfield,cutpoint);
        end

        %%%%%%%%%%%%%%%%%%%ENERGY_CALCULATION_PART%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %LOAD_WINDDATA_O1 is the first option to input the wind speed and
        %wind direction frequencies from the Winddatabase to calculate AEP.
        % In this case wind direction frequencies are listed directly by 
        % the user, while wind speed frequencies are calculated by a 
        % 2-parameters Weibull distribution.

        %% 1-基于风速和风向计算年平均发电量(风速频率则通过二参数Weibull分布来计算)
        function obj=load_winddata_o1(obj)
            obj.energy.wdata.build_fmatrix_o1();
        end

        %LOAD_WINDDATA_O2 is the second option to input the wind speed and
        %wind direction frequencies from the Winddatabase to calculate AEP.
        % In this case wind direction frequencies and wind speed 
        % frequencies are listed directly but separately by the user and
        % then combined by the software

        %% 2-基于风速和风向计算年平均发电量(分别列出风速和风向频率并结合)
        function obj=load_winddata_o2(obj)
            obj.energy.wdata.build_fmatrix_o2();
        end

        %LOAD_WINDDATA_O3 is the third option to input the wind speed and
        %wind direction frequencies from the Winddatabase to calculate AEP. 
        %This is the most preferred option due to the precision of the 
        %frequencies. In fact, in this option the frequency of each 
        %combination of wind speed and wind direction is directly reported
        %by the user. The amount of bins for velocity and direction can be 
        %chosen by the user, but we recommend 1 m/s for velocity and 5° for
        %direction

        %% 3-基于风速和风向计算年平均发电量(风速风向每个组合的频率由用户决定)
        function obj=load_winddata_o3(obj)
            obj.energy.wdata.build_fmatrix_o3();
        end
        
        %PLOT_WSPEEDS is a function that, given the frequency
        %distributions of the wind speeds, plots an histogram with the
        %relative frequency for each wind speed bin. Wind data have to be
        %loaded before calling this function.
        function obj=plot_wspeeds(obj)
            if isempty(obj.energy.wdata.frequency_matrix)
                error('Please load wind data before plotting')
            else
                obj.energy.wdata.plot_ws_distribution();
            end
        end

        %PLOT_WDIRECTIONS is a function that, given the frequency
        %distributions of the wind directions, plots a polar histogram
        %with the relative frequency for each wind direction bin. Wind 
        %data have to be loaded before calling this function.
        function obj=plot_wdirections(obj)
            if isempty(obj.energy.wdata.frequency_matrix)
                error('Please load wind data before plotting')
            else
                obj.energy.wdata.plot_wd_distribution();
            end            
        end

        %PLOT_WINDROSE is a function that, given the frequency
        %distributions of wind speeds and wind directions, plots a wind
        %rose polar histogram. Differently from the wind directions plot,
        %for each wind direction in this plot there is a subdivision in
        %wind speed frequency classes. Wind data have to be loaded before
        %calling this function.

        %% 风玫瑰图
        function obj=plot_windrose(obj)
            if isempty(obj.energy.wdata.frequency_matrix)
                error('Please load wind data before plotting')
            else
            obj.energy.wdata.plot_windrose_distribution();
            end
        end
        
        %RESET_WINDDATA is a function to empty all the data previously
        %loaded as wind speed and wind direction frequencies. Use this
        %function if it is necessary to change the input wind speed or wind
        %direction distributions or if it is necessary to change loading
        %option. 

        %% 清空风速风向数据
        function obj=reset_winddata(obj)
            obj.energy=Energy();
        end
        
        
        %CALCULATE_AEP_NOWAKE is a function that outputs the Annual Energy
        %Production in the irrealistic hypothesis of the absence of wake.
        %This calculation is useful to compute the energy efficiency of the
        %wind plant for the whole year. Basically, the "calculate_nowake"
        %function is repeated for each combination of wind speed and wind
        %direction bin (except for combinations whose frequency is equal to
        %zero) and multiplied by its respective absolute frequency during a
        %year. The results are summed together. Wind data have to be loaded 
        % before calling this function.
        
        %% 假设不存在尾流影响下的年平均发电量
        function obj=calculate_aep_nowake(obj)
            if isempty(obj.energy.wdata.frequency_matrix)
                error('Please load wind data before calculating AEP')
            else    
            obj.energy.calculate_energy_nowake(obj);
            end
        end
        

        %CALCULATE_AEP_WAKE is an important function that outputs the
        %Annual Energy Production considering the losses due to the wake
        %presence. Basically, the "calculate_wake" function is repeated for
        %each combination of wind speed and wind direction bin (except for
        %combinations whose frequency is equal to zero) and multiplied by
        %its respective absolute frequancy during a year. The results are
        %summed together. Wind data have to be loaded before calling this
        %function.

        %% 考虑尾流影响下的年平均发电量
        function obj=calculate_aep_baseline(obj)
            if isempty(obj.energy.wdata.frequency_matrix)
                error('Please load wind data before calculating AEP')
            else
           obj.energy.calculate_energy_wake(obj);
            end
        end

        %RESET_ENERGIES is a function that enables the user to delete all
        %the properties stored in the object "Energies", so that only the
        %wind data loaded are kept saved. Use this function if it is
        %necessary to delete the results of the previous energies
        %calculations without affecting the loaded wind data

        %% 删除存储在Energy中的属性
        function obj=reset_energies(obj)
            obj.energy.reset()
        end
        
        %PLOT_AEP_WSPEEDS is a function that plots a histogram with the
        %velocity bins on x-axis and with the relative contribution of each
        %velocity bin to the AEP on the y-axis. Wind data have to be
        %loaded and AEP has to be calculated before calling this function.
        function obj=plot_aep_wspeeds(obj)
            if isempty(obj.energy.aep)
                error('Please calculate AEP before plotting')
            else    
               obj.energy.plot_energy_by_speed()
            end
        end 
        
        %PLOT_AEP_WDIRECTIONS is a function that plots a histogram with the
        %direction bins on x-axis and with the relative contribution of 
        %each direction bin to the AEP on the y-axis. Wind data have to be
        %loaded and AEP has to be calculated before calling this function.
        function obj=plot_aep_wdirections(obj)
            if isempty(obj.energy.aep)
                error('Please calculate AEP before plotting')
            else    
               obj.energy.plot_energy_by_direction()
            end
        end

        %PLOT_EFF_WDIRECTIONS is a function that, given as an input a
        %velocity value, plots a polar histogram with the efficiency of
        %that wind farm for each direction. It is useful to evaluate for
        %which directions the wind farm is mostly penalized. Wind data have
        %to be loaded before calling this function.
        function obj=plot_eff_wdirections(obj,ws)
           if isempty(obj.energy.aep) || isempty(obj.energy.aep_nowake)
                error('Please calculate AEP before calculating efficiency')
           else    
              obj.energy.calculate_efficiency(obj,ws);
           end
        end

        %REPORT_ENERGIES is a function that prints in the Matlab Command
        %Window the values of the AEP considering the wake effect, AEP in
        %the hypothetical condition of no wake, the resulting efficiency
        %and the wake losses. Wind data have to be loaded and AEP have to
        %be calculated before calling this function.
        function obj=report_energies(obj)
            if isempty(obj.energy.aep) || isempty(obj.energy.aep_nowake)
                error('Please calculate AEP before printing report')
            else    
               obj.energy.report_wakeloss()
            end
        end
  
%%%%%%%%%%%%%%%%%%%%%%%%OPTIMIZATION_PART%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

        %YAW_OPTIMIZATION_GB is a function that uses the SQP algorithm
        %(gradient-based) to find the best configuration of yaw_angles to 
        %maximize power. The drawbacks of this optimization method are the
        %longer computational time if the number of turbine rises and the
        %risk to encounter local minima for some configurationsthat could 
        %not output the maximum power.

       % %% 利用基于梯度的方法来优化各风机的偏航角度
       %  function opt_yaw_angles = yaw_optimization_gb(obj)
       %      minimum_yaw_angle=obj.ya_lower; %偏航角约束
       %      maximum_yaw_angle=obj.ya_upper; %偏航角约束
       %      indexes=1:1:length(obj.layout_x); %所有风机的索引
       %      not_affecting_turbines=obj.windfield.calculate_affturbines(); %影响其他风机的风机
       %      n_turbs=length(indexes)-length(not_affecting_turbines); %影响其他风机的风机个数 
       %      if n_turbs == 0
       %         opt_yaw_angles=zeros(1,length(indexes)); %如果所有风机的尾流都不影响其他风机，那么所有风机的偏航角度都为0，返回所有风机的偏航角度
       %         return
       %      end
       %      for i=1:length(not_affecting_turbines) %从所有风机中除去不影响其他风机的索引，留下尾流对其他风机有影响的风机的索引
       %          pos=indexes==not_affecting_turbines(i);
       %          indexes(pos)=[];
       %      end
       %      opts = optimoptions('fmincon','Algorithm','sqp',...
       %          'StepTolerance',10^(-3));
       %      %x0=25*rand(1,n_turbs);
       %      x0=repelem(12,n_turbs);%12个n_turbs
       %      %x0=obj.get_yaw_angles;
       %      fun=@(x)obj.cost_function(x,indexes);
       %      A=[];
       %      b=[];
       %      Aeq=[];
       %      beq=[];
       %      lb=repelem(minimum_yaw_angle,n_turbs);%%每台风机偏航角的下界
       %      ub=repelem(maximum_yaw_angle,n_turbs);%%每台风机偏航角的上界         
       %      nonlcon=[];%没有非线性约束
       %      opt_yaw_angles_partial=fmincon(fun,x0,A,b,Aeq,beq,lb,ub,...
       %      nonlcon,opts);%求解优化问题，不含等式约束，优化变量是所有风机的偏航角度
       %      opt_yaw_angles=zeros(length(obj.layout_x),1);
       %      for i=1:length(indexes)
       %          opt_yaw_angles(indexes(i))=opt_yaw_angles_partial(i);
       %      end
       %      obj.set_yaw_angles(opt_yaw_angles);
       %      %opt_yaw_angles=round(opt_yaw_angles);
       %      %obj.set_yaw_angles(opt_yaw_angles);
       %      obj.calculate_wake();
       %  end

        %% 利用基于梯度的方法来优化各风机的轴向诱导因子
        function opt_axial_induction_factor = aif_optimization_gb(obj)
            minimum_aif=obj.aif_lower; %偏航角约束
            maximum_aif=obj.aif_upper; %偏航角约束
            indexes=1:1:length(obj.layout_x); %所有风机的索引
            % not_affecting_turbines=obj.windfield.calculate_affturbines(); %影响其他风机的风机
            % n_turbs=length(indexes)-length(not_affecting_turbines); %影响其他风机的风机个数 
            % if n_turbs == 0
            %    opt_axial_induction_factor=zeros(1,length(indexes)); %如果所有风机的尾流都不影响其他风机，那么所有风机的偏航角度都为0，返回所有风机的偏航角度
            %    return
            % end
            % for i=1:length(not_affecting_turbines) %从所有风机中除去不影响其他风机的索引，留下尾流对其他风机有影响的风机的索引
            %     pos=indexes==not_affecting_turbines(i);
            %     indexes(pos)=[];
            % end
            opts = optimoptions('fmincon','Algorithm','sqp',...
                'StepTolerance',10^(-3));
            %x0=25*rand(1,n_turbs);
            x0=repelem(1/3,67);%12个n_turbs
            
            %x0=obj.get_yaw_angles;
            fun=@(x)obj.cost_function(x,indexes);
            A=[];
            b=[];
            Aeq=[];
            beq=[];
            lb=repelem(minimum_aif,67);%%每台风机偏航角的下界
            ub=repelem(maximum_aif,67);%%每台风机偏航角的上界         
            nonlcon=[];%没有非线性约束
            opt_axial_induction_factor_partial=fmincon(fun,x0,A,b,Aeq,beq,lb,ub,...
            nonlcon,opts);%求解优化问题，不含等式约束，优化变量是所有风机的偏航角度
            opt_axial_induction_factor=zeros(length(obj.layout_x),1);
            for i=1:length(indexes)
                opt_axial_induction_factor(indexes(i))=opt_axial_induction_factor_partial(i);
            end
            obj.set_aif(opt_axial_induction_factor);
            %opt_axial_induction_factor=round(opt_axial_induction_factor);
            %obj.set_yaw_angles(opt_axial_induction_factor);
            obj.calculate_wake();
        end





        %aif_optimization_ga is a function that uses the genetic algorithm
        %to find the best configuration of yaw angles to maximize power.
        %The integer constraint is done mainly to reduce the computational
        %time and due to the fact that decimal yaw angles are not of
        %practical interest. Max generation number is set to 200, while
        %Population size to 10; this leads to acceptable computational time
        %and satisfactory results. However, since GA is a stochastic
        %algorithm, it may happen that at the end of a simulation some
        %yaw angles have some errors and power is not completely optimized.
        %A higher number of populations would reduce even further this
        %possibility but would lead to unacceptable computational time.

        %% 
        function opt_axial_induction_factor = aif_optimization_ga(obj)
            minimum_aif=obj.aif_lower;
            maximum_aif=obj.aif_upper;
            indexes=1:1:length(obj.layout_x);
            not_affecting_turbines=obj.windfield.calculate_affturbines();
            n_turbs=length(indexes)-length(not_affecting_turbines);
            if n_turbs == 0
               opt_axial_induction_factor=zeros(1,length(indexes));
               return
            end
            for i=1:length(not_affecting_turbines)
                pos=indexes==not_affecting_turbines(i);
                indexes(pos)=[];
            end
            rng default
            fun=@(x)obj.cost_function(x,indexes);
            nvars=n_turbs;
            Aineq=[];
            Bineq=[];
            Aeq=[];
            beq=[];
            lb=repelem(minimum_aif,n_turbs);
            ub=repelem(maximum_aif,n_turbs);            
            nonlcon=[];
            intcon=linspace(1,nvars,nvars);
            genmax=200+length(obj.layout_x);
            opts = optimoptions('ga','MaxGenerations',genmax,...
            'PopulationSize',12,'PlotFcn', @gaplotbestf);
            opt_axial_induction_factor_partial=ga(fun,nvars,Aineq,Bineq,Aeq,beq,lb,ub,...
            nonlcon,intcon,opts);
            opt_axial_induction_factor=zeros(length(obj.layout_x),1);
            for i=1:length(indexes)
                opt_axial_induction_factor(indexes(i))=opt_axial_induction_factor_partial(i);
            end
            obj.set_aif(opt_axial_induction_factor);
            obj.calculate_wake();
        end

        %YAW_OPTIMIZATION_MIGA is a function that uses Mixed Integer GA
        %Optimization to determine the optimal configuration of yaw angles.
        %Here, the values that yaw angles could assume to optimize total
        %farm power are restricted to some equally spaced values between
        %the admitted boundaries. This is an interesting option since in
        %reality is not required a precision of 1° of yaw angles due to
        %measurement uncertainties and the negligible effect of finely
        %adjusting the yaw angles on total power. This option lets to
        %decrease the parameter Max Generation and as a result reduce the
        %compuational time. However, since GA is a stochastic algorithm, it
        %may happen that at the end of a simulation some yaw angles have
        %some errors and power is not completely optimized.
        function opt_axial_induction_factor = aif_optimization_miga(obj)
            n_turbs=length(obj.layout_x);  
            low=obj.aif_lower;
            upp=obj.aif_upper;
            fun=@(x)obj.fast_cost_function(x,low,upp);
            nvars=n_turbs;
            Aineq=[];
            Bineq=[];
            Aeq=[];
            beq=[];
            lb=repelem(1,n_turbs);
            ub=repelem(obj.ya_options,n_turbs);            
            nonlcon=[];
            intcon=linspace(1,nvars,nvars);
            opts = optimoptions('ga','MaxGenerations',130,...
                'PopulationSize',10,'PlotFcn', @gaplotbestf);
            opt_axial_induction_factor=ga(fun,nvars,Aineq,Bineq,Aeq,beq,lb,ub,...
                nonlcon,intcon,opts);
            opt_axial_induction_factor=low+(upp-low)*(opt_axial_induction_factor-1)...
                /(obj.ya_options-1);
            obj.set_aif(opt_axial_induction_factor);
        end

        %YAW_OPTIMIZATION_SQ is a function that optimizes the yaw 
        %angles to maximize total farm power. This algorithm is
        %gradient-free and assumes that the optimal angle of a turbine is
        %not correlated with the yaw angles of the downstream turbines.
        %For each turbine i sequentially, all turbine
        %options are simulated and the farm power for each option is stored
        %in a matrix. The yaw angle that maximizes the farm power is chosen
        %and fixed for that turbine. The process is repeated for all other
        %turbines. The number ofoptions that each turbine can assume can be
        %chosen by the user. This algorithm gives a speed computational
        %advantage with respect to other algorithms since it is basically a
        %sequence of calculate_wake(). Furthermore, from the tests made
        %with several layouts, it is also the most precise in finding the
        %optimal yaw angles.
        function opt_axial_induction_factor=aif_optimization_sq(obj)
            n_turbs=length(obj.layout_x);  
            aif_t=zeros(n_turbs,1);
            winning_power=0;
            options=linspace(obj.aif_lower,obj.aif_upper,obj.ya_options);
            p=zeros(length(options),1);
            sorted_indexes=obj.windfield.get_ordered_turbines();
            not_affecting_turbines=obj.windfield.calculate_affturbines();
            if length(sorted_indexes)==length(not_affecting_turbines)
               opt_axial_induction_factor=zeros(1,length(sorted_indexes));
               return
            end
            for i=1:length(not_affecting_turbines)
                pos=sorted_indexes==not_affecting_turbines(i);
                sorted_indexes(pos)=[];
            end
            if any(obj.turbine_status==0)
                fau=find(~obj.turbine_status);
                for i=1:length(fau)
                    idx=sorted_indexes==fau(i);
                    sorted_indexes(idx)=[];
                end
            end
            for j=1:length(sorted_indexes)
                for k=1:length(options)
                    aif_t(sorted_indexes(j))=options(k);
                    obj.set_aif(aif_t);
                    obj.calculate_wake();
                    p(k)=obj.get_farm_power();
                end
                if max(p)>winning_power
                    [~,fau]=max(p);
                    aif_t(sorted_indexes(j))=options(fau);
                end
            end
            opt_axial_induction_factor=aif_t';
            obj.set_aif(opt_axial_induction_factor);
            obj.calculate_wake();
        end

%ANN_optimizer, loading file is necessary
%         function [opt_axial_induction_factor,final_power]=nn_optimizer(obj)
%             load Mdl.mat;
%             a=obj.windfield.wind_speed;
%             b=obj.windfield.wind_direction;
%             n_turbs=length(obj.layout_x);  
%             yaw_angles=zeros(n_turbs,1);
%             winning_power=0;
%             options=linspace(obj.aif_lower,obj.ya_upper,obj.ya_options);
%             p_matrix=zeros(length(options),1);
%             sorted_indexes=obj.windfield.get_ordered_turbines();
%             not_affecting_turbines=obj.windfield.calculate_affturbines();
%             for i=1:length(not_affecting_turbines)
%                 pos=sorted_indexes==not_affecting_turbines(i);
%                 sorted_indexes(pos)=[];
%             end
%             for j=1:length(sorted_indexes)
%                 for k=1:length(options)
%                     yaw_angles(sorted_indexes(j))=options(k);
%                     c=yaw_angles;
%                     p_matrix(k)=predict(Mdl,[a b c']);
%                 end
%                 if max(p_matrix)>winning_power
%                     [winning_power,idx]=max(p_matrix);
%                     yaw_angles(sorted_indexes(j))=options(idx);
%                 end
%             end
%             opt_axial_induction_factor=yaw_angles';
%             obj.set_aif(opt_axial_induction_factor);
%             d=opt_axial_induction_factor;
%             final_power=predict(Mdl,[a b d]);
%         end
        
        %CALCULATE_AEP_OPTIMIZED is a function that outputs the AEP
        %optimized due to yaw steering for a given relative frequancy
        %distribution of wind speeds and wind directions. Wind data have to
        %be loaded before calling this function.
        function obj=calculate_aep_optimized(obj)
            if isempty(obj.energy.wdata.frequency_matrix)
                error('Please load wind data before calculating AEP opt')
            else
                obj.energy.energy_opt(obj)
            end
        end

        %REPORT_ENERGIES_OPT is a function that prints in the Matlab 
        %Command Window the values of the AEP considering the wake effect,
        %AEP optimized, the increased percentage in AEP, and the additional
        %revenue. AEP optimized has to be calculated before calling this
        %function.
        function obj=report_energies_opt(obj)
            if isempty(obj.energy.aep_opt)
                error('Please calculate opt AEP before printing report')
            else    
               obj.energy.report_yaw_optimization()
            end
        end

        %PLOT_GAIN_WDIRECTIONS is a function that, given a wind_speed as an
        %input, returns a polar histogram with the gain in power for each
        %wind direction after yaw optimization. This function is effective
        %to investigate which directions are more affected by the yaw
        %steering than others. Optimized AEP has to be calculated before
        %calling this function.
        function obj=plot_gain_wdirections(obj,ws)
           if isempty(obj.energy.aep_opt)
                error('Please calculate opt AEP before printing report')
           else    
              obj.energy.calculate_relative_power_gain(ws);
           end
        end
    end

    methods (Access=private)
        function power=cost_function(obj,aif_t,aff_turbines)
            vector=zeros(length(obj.layout_x));
            for i=1:length(aff_turbines)
                vector(aff_turbines(i))=aif_t(i);
            end
            obj.set_aif(vector);
            obj.calculate_wake();
            power=-obj.get_farm_objective();
        end

        function power=fast_cost_function(obj,aif_t,low,upp)
            all_yaws=linspace(low,upp,obj.ya_options);
            aif_t=all_yaws(aif_t);
            obj.set_aif(aif_t);
            obj.calculate_wake();
            power=-obj.get_farm_objective();
        end    
    end
end