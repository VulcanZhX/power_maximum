classdef Wake < handle

    properties
        velocity_model=readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B3:B3');   %从excel文件'inputs.xlsx'中读取的风速模型，尾流速度计算的方法
        deflection_model=readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B4:B4'); %从excel文件'inputs.xlsx'中读取的尾流偏转模型，尾流偏转计算的方法
        combination_model=readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B5:B5');%从excel文件'inputs.xlsx'中读取的尾流组合模型，选择如何组合多个尾流速度亏损的方法
        turbulence_combination_model=readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B7:B7');%从excel文件'inputs.xlsx'中读取的尾流组合模型，选择如何组合多个尾流速度亏损的方法
        turbulence_model=readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B5:B5'); %从excel文件'inputs.xlsx'中读取的湍流模型，选择湍流计算方法
        wd=Wakedeflection  %%尾流偏转模型
        wv=Wakevelocity    %%尾流风速模型
        wc=Wakecombination %%尾流叠加模型
        wt=Waketurbulence  %%尾流湍流模型
        ti_initial=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B9:B9'));%湍流强度参数
        ti_constant=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B10:B10'));%湍流强度参数
        ti_ai=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B11:B11'));%湍流强度参数
        ti_downstream=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B12:B12'));%湍流强度参数
        turbulence_choice=readcell('inputs_all_fields.xlsx','Sheet','Models','Range','F2:F2');%是否考虑湍流
    end

    methods
        function [deflection]=deflection_function(obj,x_locations,y_locations,turbine,coord,wind_field,u_initial)%风机尾流偏转计算，
            if strcmp(obj.deflection_model,'Huadian')
                [deflection]=obj.wd.Huadian_deflect_function_yaw(x_locations,y_locations,turbine,coord,wind_field,u_initial);
            elseif strcmp(obj.deflection_model,'Gauss')
                [deflection]=obj.wd.Gauss_function(x_locations,y_locations,turbine,coord,wind_field,u_initial);
            elseif strcmp(obj.deflection_model,'Jimenez')
                [deflection]=obj.wd.Jimenez_function(x_locations,y_locations,turbine,coord,wind_field,u_initial);
            else
                error('Input deflection model not valid');
            end
        end

        function [deficit1,deficit2,deficit3]=velocity_function(obj,x_locations,y_locations,z_locations,turbine,coord,deflection,wind_field,u_initial)
            if strcmp(obj.velocity_model,'Jensen')
                [deficit1,deficit2,deficit3]=obj.wv.Jensen_function(x_locations,y_locations,z_locations,turbine,coord,deflection,wind_field,u_initial);
            elseif strcmp(obj.velocity_model,'Multizone')
                [deficit1,deficit2,deficit3]=obj.wv.Multizone_function(x_locations,y_locations,z_locations,turbine,coord,deflection,wind_field,u_initial);
            elseif strcmp(obj.velocity_model,'Gauss')
                [deficit1,deficit2,deficit3]=obj.wv.Gauss_function(x_locations,y_locations,z_locations,turbine,coord,deflection,wind_field,u_initial);
            elseif strcmp(obj.velocity_model,'Huadian')
                [deficit1,deficit2,deficit3]=obj.wv.Huadian_velocity_function_yaw(x_locations,y_locations,z_locations,turbine,coord,deflection,wind_field,u_initial);
            else
                error('Input velocity model not valid');
            end
        end

        function deficit_combination=combination_function(obj,field,wake)
            if strcmp(obj.combination_model,'Linear')
                deficit_combination=obj.wc.linear_function(field,wake);
            elseif strcmp(obj.combination_model,'Sos')
                deficit_combination=obj.wc.sumofsquares_function(field,wake);
            else
                error('Input combination model not valid');
            end
        end

        function deficit_combination_turbulence=turbulence_combination_function(obj,field,wake)
            if strcmp(obj.combination_model,'Linear')
                deficit_combination_turbulence=obj.wc.linear_function(field,wake);
            elseif strcmp(obj.combination_model,'Sos')
                deficit_combination_turbulence=obj.wc.sumofsquares_function(field,wake);
            else
                error('Input combination model not valid');
            end
        end

        
        
        function wake_turbulence=turbulence_function(obj,x_locations,y_locations,z_locations,turbine,coord,deflection,windfield,turbulence_init)
           if strcmp(obj.turbulence_model,'Huadian') 
               wake_turbulence=obj.wt.Huadian_turbulence_function_yaw(x_locations,y_locations,z_locations,turbine,coord,deflection,windfield,turbulence_init);
           else
           end

        end
    end
end