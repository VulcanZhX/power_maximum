classdef  Waketurbulence< handle

    properties
        Sc_t=0.5;
        S_1 = 0.043;
        sigma_e = 0.18;
        wake_width
        %wake_turbulence


    end

    methods

            %% 现场尾流湍流计算模型(华电葛老师提供)
                function [wake_turbulence]=Huadian_turbulence_function(obj,x_locations,y_locations,z_locations,turbine,coord,deflection,windfield,turbulence_init)
        
                %% 沿流向，各待计算点的距离
                x_d=x_locations-coord(1);
                I_d=zeros(size(x_d)); % 各待计算点湍流强度初始化
                mask = x_d >= 0; % 只计算下游点
        
                %沿展向，各待计算点的距离
                y_d=y_locations-coord(2);
                y_d(abs(y_d) < 0.0001) = 0.0001; 
        
                %沿垂向，各待计算点的距离
                z_d=z_locations-coord(3);
        
                %% 指数项
                dr=sqrt((y_d - deflection).^2 + z_d.^2);
                % dr
        

                %% 具体湍流强度计算

                % 用到的参数
                D=turbine.rotor_diameter;
                gamma=turbine.yaw_angle;
                Ct=turbine.Ct;
                Ia=turbine.turbulence;

                % 推力系数接近0，附加湍流强度为0
                
                d = 2.3*(Ct*cosd(gamma))^(-1.2);
                e = Ia^0.1;
                q = 0.7*(Ct*cosd(gamma))^(-3.2)*Ia^(-0.45).*(1 + x_d./D).^(-2);
                dImax = 1./(d+ e.*x_d./D+q);
        
                
                %% 计算尾流宽度  
                obj.calculate_wake_width(x_locations,turbine,coord,windfield);
                WD = obj.wake_width;
                

                %% 标准尺度换算
                sigma_y=WD*D;
                
                sigma_tm=sigma_y./sqrt(2*log(2));
                R_half=sqrt(2*log(2)).*sigma_y;
                % R_half
        
                % 横向分布
                shape_tm=zeros(size(dr));
                k1=zeros(size(dr));
        
                % Case1: dr<R_half
                cond=dr<R_half;
                shape_tm(cond)=1-0.15*(1 + cosd(pi*dr(cond)./R_half(cond)));
                k1(cond) = sind((pi / 2) * dr(cond) ./ R_half(cond));
        
                % Case2: dr>=R_half
                shape_tm(~cond)=exp(-((dr(~cond)-R_half(~cond)).^2)./(2 * sigma_tm(~cond).^2));
                k1(~cond)=1;
                alpha_tm=atan2d(abs(z_d), abs(y_d));
        
                % 垂向修正
                delta_tm=zeros(size(z_d));
                exp_part=exp(-((dr-R_half).^2)./(2*sigma_tm.^2));
                delta_tm(z_d>=0)=0.23*Ia.*sind(alpha_tm(z_d>=0)).*(k1(z_d>=0).* exp_part(z_d>=0));
                delta_tm(z_d<0)=-1.23*Ia.*sind(alpha_tm(z_d<0)).*(k1(z_d<0).*exp_part(z_d<0));
        
                % 结果
                I_d(mask)=dImax(mask).*shape_tm(mask)+delta_tm(mask);
                wake_turbulence=I_d.*turbulence_init;

                end
            

                function calculate_wake_width(obj,x_locations,turbine,coord,windfield)
                %% 沿流向，各待计算点的距离
                x_d=x_locations-coord(1);
                
        
               
        
                
        
               
        

                %% 具体湍流强度计算

                % 用到的参数
                D=turbine.rotor_diameter;
                Ct=turbine.Ct;
                Ia=turbine.turbulence;                           % 自然风下的湍流强度

        
                
                beta_1=0.5 * (1+sqrt(1-Ct))/sqrt(1-Ct);
                epsilon=0.2 *sqrt(beta_1);
                k_w = 0.38 * Ia + 0.004; %%尾流膨胀系数


                %% 计算近尾流
                if abs(windfield.wind_speed - turbine.average_velocity) < 0.005     
                x_0=D;
                else
                x_0=0.5*D;
                end
                x_NW =x_0+obj.sigma_e*D*(1+sqrt(1-Ct))/(2*(sqrt(obj.Sc_t)*(0.63*Ia)+obj.S_1*(1-sqrt(1-Ct))));
                
                
                %% 计算尾流宽度
                %1) 计算 delta_Umax
                delta_Umax = zeros(size(x_d));
                idx1 = x_d<D;
                idx2 = ~idx1;
                
                delta_Umax(idx1) = 1.75 *(D/x_NW +0.5).^(-1.37).*(1-sqrt(1-Ct));
                delta_Umax(idx2) = 1.75 *(x_d(idx2)/x_NW + 0.5).^(-1.37).*(1-sqrt(1-Ct));
                
                % 2) 计算K
                K = zeros(size(x_d));
                idx1 = x_d<8*D;
                idx2 = (x_d>=8*D)&(x_d<=12*D);
                idx3 = x_d>12*D;
                K(idx1) = 1;
                K(idx2) = 1 - (x_d(idx2)-8*D)/(4*D);
                K(idx3) = 0;  % 这行可省略，因为初始化为 0
    
                % 3)计算尾流宽度
                 obj.wake_width = K.*(k_w.*x_d/D+epsilon)+(1-K).*sqrt(Ct./(8*(1-(1-delta_Umax).^2)));
                end

                
                %% 含偏航的尾流宽度计算
                function calculate_wake_width_yaw(obj,x_locations,turbine,coord,~)
                %% 沿流向，各待计算点的距离
                x_d=x_locations-coord(1);
                rotor_D=turbine.rotor_diameter;
                obj.wake_width =1 + 0.0834 * log(1 + exp((x_d - rotor_D) ./ rotor_D * 2));
                end

                %% 含偏航的湍流计算
              function [wake_turbulence]=Huadian_turbulence_function_yaw(obj,x_locations,y_locations,z_locations,turbine,coord,deflection,windfield,turbulence_init)
        
                %% 沿流向，各待计算点的距离
                x_d=x_locations-coord(1);
                I_d=zeros(size(x_d)); % 各待计算点湍流强度初始化
                mask = x_d >= 0; % 只计算下游点
        
                %沿展向，各待计算点的距离
                y_d=y_locations-coord(2);
                y_d(abs(y_d) < 0.0001) = 0.0001; 
        
                %沿垂向，各待计算点的距离
                z_d=z_locations-coord(3);
        
                %% 指数项
                dr=sqrt((y_d - deflection).^2 + z_d.^2);
                
        

                %% 具体湍流强度计算

                % 用到的参数
                D=turbine.rotor_diameter;
                gamma=turbine.yaw_angle;
                Ct=turbine.Ct;
                Ia=turbine.turbulence;

                % 推力系数接近0，附加湍流强度为0
                
                d = 2.3*(Ct*cosd(gamma))^(-1.2);
                e = Ia^0.1;
                q = 0.7*(Ct*cosd(gamma))^(-3.2)*Ia^(-0.45).*(1 + x_d./D).^(-2);
                dImax = 1./(d+ e.*x_d./D+q);
        
                
                %% 计算尾流宽度  
                obj.calculate_wake_width_yaw(x_locations,turbine,coord,windfield);
                WD = obj.wake_width;
                

                %% 标准尺度换算
                sigma_0=0.25*D;
                sigma_y=sigma_0.*WD;
                
                sigma_tm=sigma_y./sqrt(2*log(2));
                R_half=sqrt(2*log(2)).*sigma_y;
                % R_half
        
                % 横向分布
                shape_tm=zeros(size(dr));
                k1=zeros(size(dr));
        
                % Case1: dr<R_half
                cond=dr<R_half;
                shape_tm(cond)=1-0.15*(1 + cosd(180*dr(cond)./R_half(cond)));
                k1(cond) = sind(90 * dr(cond) ./ R_half(cond));
        
                % Case2: dr>=R_half
                shape_tm(~cond)=exp(-((dr(~cond)-R_half(~cond)).^2)./(2 * sigma_tm(~cond).^2));
                k1(~cond)=1;
                alpha_tm=atan2d(abs(z_d), abs(y_d));
        
                % 垂向修正
                delta_tm=zeros(size(z_d));
                exp_part=exp(-((dr-R_half).^2)./(2*sigma_tm.^2));
                delta_tm(z_d>=0)=0.23*Ia.*sind(alpha_tm(z_d>=0)).*(k1(z_d>=0).* exp_part(z_d>=0));
                delta_tm(z_d<0)=-1.23*Ia.*sind(alpha_tm(z_d<0)).*(k1(z_d<0).*exp_part(z_d<0));
        
                % 结果
                I_d(mask)=dImax(mask).*shape_tm(mask)+delta_tm(mask);
                wake_turbulence=I_d.*turbulence_init;

              end



                  
     
   
    end
end