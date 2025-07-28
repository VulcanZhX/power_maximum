classdef Turbinechart < handle

    properties
        layout_x                                                                     %风机布局图的x的坐标数组
        layout_y                                                                     %风机布局图的y的坐标数组
        turbines
    end

    properties (Dependent)
        layout_z                                                                     %风机布局图的z坐标数组
        coordinates                                                                  %风机的坐标数组  
        indexes                                                                      %风机的索引数组
        turbines_cell_array                                                          %风机信息的单元格数组
    end

    methods

        function obj=Turbinechart(layout_x,layout_y,turbines)
            obj.layout_x=layout_x;                %存储所有风机的x坐标
            obj.layout_y=layout_y;                %存储所有风机的y坐标
            obj.turbines=turbines;                %风机信息的单元格数组
        end

        
        %% 获取风机布置图中每台风机的z坐标
        function layout_z=get.layout_z(obj)
            layout_z=zeros(length(obj.layout_x),1);
            for i=1:length(obj.layout_x)
                layout_z(i,1)=obj.turbines{i,1}.hub_height;   %trubines{i,1}第i台风机的轮毂高度
            end
        end

        %% 获取风机布局图中每台风机的三维坐标信息
        function coordinates=get.coordinates(obj)
            coordinates=[obj.layout_x,obj.layout_y,obj.layout_z];
        end
        
        %% 当前风场布局下，每台风机的顺序索引
        function indexes=get.indexes(obj)
            indexes=(1:1:length(obj.layout_x))';
        end 

        %% cell数组包含所有风机的坐标，机组信息以及对应的索引
        function turbines_cell_array=get.turbines_cell_array(obj)  %cell用于存储整个wind farm中所有风机的坐标、机组以及索引信息
            turbines_cell_array=cell(length(obj.layout_x),3);
            for i=1:length(obj.layout_x)
            turbines_cell_array{i,1}=obj.coordinates(i,:);         %turbines_cell_array表示一行一风机，所有cell行的第一列存储坐标信息
            turbines_cell_array{i,2}=obj.turbines{i,1};            %所有cell行的第二列存储风机对象信息，风机对象类
            turbines_cell_array{i,3}=obj.indexes(i,1);             %所有cell行的第三列存储风机索引信息
            end
        end
    
        %% 风向变化后，风机对应的旋转坐标
        function rotated_turbines_cell_array=rotated(obj,center,angle)
            cell_array=obj.turbines_cell_array;
            rotated_turbines_cell_array=cell_array;
            for i=1:length(obj.layout_x)
                point=[cell_array{i,1}(1),cell_array{i,1}(2),cell_array{i,1}(3)];   %提取当前风机的坐标点
                [r]=rotation(point,center,angle);                                   %风向改变后，风机的旋转坐标
                rotated_turbines_cell_array{i,1}(1)=r(1);                           %风向改变后，风机旋转对应的x坐标
                rotated_turbines_cell_array{i,1}(2)=r(2);                           %风向改变后，风机旋转对应的y坐标
                rotated_turbines_cell_array{i,1}(3)=r(3);                           %风向改变后，风机旋转对应的z坐标
            end
        end
    end

    methods(Static)
        
        %% 对风机的cell数组进行重新排序，风机的顺序索引是不变的
        function sorted_turbines_cell_array=sortinx(input_cell)
            [l,~]=size(input_cell);           %cell数组的行大小
            indexes=zeros(l,1);
            for i=1:l
                indexes(i)=input_cell{i,1}(1);%输入数组所有行第一列数组第一个元素，index为各风机在x方向上的坐标
            end
            [~,indexes_sorted]=sort(indexes); %依据x的坐标最风机由小到大进行排序，重新排序后原风机的索引不变
            sorted_turbines_cell_array=cell(size(input_cell));                  %创建sorted_turbines_cell_array，和原来的input_cell大小一致
            for i=1:l
                sorted_turbines_cell_array{i,1}=input_cell{indexes_sorted(i),1};%重新排序后的风机cell数组，按x坐标的大小顺序排序，
                sorted_turbines_cell_array{i,2}=input_cell{indexes_sorted(i),2};
                sorted_turbines_cell_array{i,3}=input_cell{indexes_sorted(i),3};                
            end
        end
    end
end