classdef Wakedeflection < handle

    properties
        kd=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B33:B33'));%Jimenez湍流模型参数
        ad_j=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B34:B34'));%Jimenez 模型中的额外参数
        bd_j=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B35:B35'));%Gauss 模型中使用的参数。
        ka=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B25:B25'));
        kb=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B26:B26'));
        alpha=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B27:B27'));
        beta=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B28:B28'));
        ad_g=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B29:B29'));
        bd_g=cell2mat(readcell('inputs_all_fields.xlsx','Sheet','Models','Range','B30:B30'));
    end

    methods
        function [deflection_jimenez]=Jimenez_function(obj,x_locations,~,turbine,coord,~,~)
            rotor_diameter=turbine.rotor_diameter;%风机转子直径
            yaw_angle=-turbine.yaw_angle;
            Ct0=(cosd(yaw_angle))^turbine.pT;
            aif=turbine.axial_induction_factor;%风机推力系数
            alpha_0=(cosd(yaw_angle))^2*sind(yaw_angle)*Ct0*aif/2;
            x_new=x_locations-coord(1);
            %alpha=rad2deg(alpha_0*(1/(1+(obj.kd*x_new/turbine.rotor_radius))).^2);
            dx=alpha_0*(15*(2*obj.kd*x_new/rotor_diameter+1).^4+alpha_0^2)...
            ./(30*obj.kd/rotor_diameter*(2*obj.kd*x_new/rotor_diameter+1).^5)-...
            alpha_0*rotor_diameter*(15+alpha_0^2)/30/obj.kd;
            deflection_jimenez=-(dx+obj.ad_j+obj.bd_j*x_new);
        end
        

        %% 现场尾流偏转模型(简化版)
        function [deflection_huadian_deflect]=Huadian_deflect_function(~,x_locations,~,turbine,coord,~,~)
            rotor_diameter=turbine.rotor_diameter;%风机转子直径
            yaw_angle=-turbine.yaw_angle;
            % yaw_angle_rad = deg2rad(yaw_angle);
            x_new=x_locations-coord(1);
            deflection_huadian_deflect=(cosd(yaw_angle))^2*sind(yaw_angle)*turbine.Ct*5*(1-1/(1+0.1*(x_new/rotor_diameter)))*rotor_diameter;

            % Ct0=(cosd(yaw_angle))^turbine.pT;
            % aif=turbine.axial_induction_factor;%风机推力系数
            % alpha_0=(cosd(yaw_angle))^2*sind(yaw_angle)*Ct0*aif/2;
            % 
            % %alpha=rad2deg(alpha_0*(1/(1+(obj.kd*x_new/turbine.rotor_radius))).^2);
            % dx=alpha_0*(15*(2*obj.kd*x_new/rotor_diameter+1).^4+alpha_0^2)...
            % ./(30*obj.kd/rotor_diameter*(2*obj.kd*x_new/rotor_diameter+1).^5)-...
            % alpha_0*rotor_diameter*(15+alpha_0^2)/30/obj.kd;
            % deflection_jimenez=-(dx+obj.ad_j+obj.bd_j*x_new);
        end

        %% 现场尾流偏转模型(偏航版)
        function [deflection_huadian_deflect]=Huadian_deflect_function_yaw(~,x_locations,~,turbine,coord,~,~)
            %% 1.计算尾流宽度
            rotor_D=turbine.rotor_diameter;
            x_d=x_locations-coord(1);
            y_c=zeros(size(x_d));
            for i = 1:length(x_d)
                dx_val = x_d(i);
                if dx_val > 0
                    delta_x_steps = 0:1:(dx_val - 1);  % 步长为10
                    if ~isempty(delta_x_steps)
                        d_w_steps = 1 + 0.0834 * log(1 + exp((delta_x_steps - rotor_D) / rotor_D * 2));
                        integrand = -turbine.Ct * cosd(turbine.yaw_angle)^2 * sind(turbine.yaw_angle) ./ (8 * d_w_steps.^2) ...
                                    .* (1 + erf(sqrt(2) * delta_x_steps / rotor_D));
                        y_c(i) = sum(integrand) * 2;  % 矩形法积分，步长10 * 2 = 20
                    end
                else
                    y_c(i) = 0;
                end
            end
            deflection_huadian_deflect=y_c;
        end





        function [deflection_gauss]=Gauss_function(obj,x_locations,y_locations,turbine,coord,wind_field,u_initial)
            yaw_angle=-turbine.yaw_angle;
            Ct0=(cosd(yaw_angle))^turbine.pT;
            aif=turbine.axial_induction_factor;%风机推力系数
            rotor_diameter=turbine.rotor_diameter;
            turbulence=turbine.turbulence;
            wind_veer=wind_field.wind_veer;
            u_inf=u_initial;
            ur=u_inf*Ct0*aif*cosd(yaw_angle)/2/(1-sqrt(1-Ct0*aif*cosd(yaw_angle)));
            u0=u_inf*sqrt(1-Ct0*aif);
            C0=1-u0/wind_field.wind_speed;
            M0=C0.*(2-C0);
            E0=C0.^2-3*exp(1/12)*C0+3*exp(1/3);
            x0=rotor_diameter*cosd(yaw_angle)*(1+sqrt(1-Ct0*aif))/sqrt(2)/(4*obj.alpha*turbulence+2*obj.beta*(1-sqrt(1-Ct0*aif)))+coord(1);
            ky=obj.ka*turbulence+obj.kb;
            kz=obj.ka*turbulence+obj.kb;
            sigma_z0=0.5*rotor_diameter*sqrt(ur./(u_inf+u0));
            sigma_y0=sigma_z0*cosd(yaw_angle)*cosd(wind_veer);
            theta_c0=rad2deg(0.3*deg2rad(yaw_angle)*(1-sqrt(1-Ct0*aif*cosd(yaw_angle)))/cosd(yaw_angle));
            delta0=(x0-coord(1))*tand(theta_c0);
            yR=y_locations-coord(2);
            xR=yR*tand(yaw_angle)+coord(1);
            delta_nearwake=((x_locations-xR)./(x0-xR)*delta0)+obj.ad_g+obj.bd_g*(x_locations-coord(1));
            for i=1:numel(delta_nearwake)
                if x_locations(i)<xR(i) || x_locations(i)>x0
                delta_nearwake(i)=0;
                end
            end
            sigma_y=ky*(x_locations-x0)+sigma_y0;
            sigma_z=kz*(x_locations-x0)+sigma_z0;
            for i=1:numel(x_locations)
                if x_locations(i)<x0
                   sigma_y(i)=sigma_y0(i);
                   sigma_z(i)=sigma_z0(i);
                end
            end
            ln_numerator=(1.6+sqrt(M0)).*(1.6*sqrt(sigma_y.*sigma_z./sigma_y0./sigma_z0)-sqrt(M0));
            ln_denominator=(1.6-sqrt(M0)).*(1.6*sqrt(sigma_y.*sigma_z./sigma_y0./sigma_z0)+sqrt(M0));
            delta_farwake=delta0+deg2rad(theta_c0)*E0/5.2.*sqrt(sigma_y0.*sigma_z0/ky/kz./M0).*log(ln_numerator./ln_denominator)+obj.ad_g+obj.bd_g*(x_locations-coord(1));
            for i=1:numel(delta_farwake)
                if x_locations(i)<=x0
                delta_farwake(i)=0;
                end
            end
            deflection_gauss=delta_nearwake+delta_farwake;
        end
    end
end