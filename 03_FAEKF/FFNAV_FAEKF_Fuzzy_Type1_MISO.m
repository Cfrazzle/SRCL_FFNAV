function [lambda] = FFNAV_EKF_Fuzzy_Type1_MISO(x, FUZZY_props)
% FFNAV Fuzzy Controller Type 1 - MISO ====================================
% Description: This function completes the fuzzification, implication, and
% defuzzification step of the FLS. Given a two normalized inputs, the 
% function returns the corresponding normalized output of the fuzzy system.
%
% Inputs:
%   x            - Normalized inputs to the fuzzy system (two inputs) 
%   FUZZY_props. - Structure containing properties of the FLS
%
% Outputs:
%   lambda       - Normalized output of the fuzzy system
%
% Other Function Calls:
%   None
%
% Created by:  Cory Fraser - JAN 08, 2018
% Latest Edit: Cory Fraser - JUL 04, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% Initialize Parameters
x1 = x(1); 
x2 = x(2); 

rules         = FUZZY_props.rules_MISO;
n_pts         = FUZZY_props.n_pts;
X_pts         = FUZZY_props.X_pts;
Y_pts         = FUZZY_props.Y_pts;
out_mf        = FUZZY_props.out_mf_MISO;
implied_mf    = FUZZY_props.implied_mf_MISO;
aggregate_mf  = FUZZY_props.aggregate_mf;
c_IN           = FUZZY_props.c_IN;
cMISO         = FUZZY_props.c_MISO;
implied_flag  = FUZZY_props.implied_flag;
n_LingValsIn  = FUZZY_props.n_LingValsIn;
n_LingValsOut = FUZZY_props.n_LingValsOut_MISO;
sigma_IN      = FUZZY_props.sigma_IN; 
sigmoid_a       = FUZZY_props.sigmoid_a     ;
sigmoid_c       = FUZZY_props.sigmoid_c;

% =========================================================================
%% Fuzzification (Eq. 4.45, 4.46)

u_x1 = zeros(1,n_LingValsIn);

u_x1(1) = 1./(1+exp(sigmoid_a*(x1+sigmoid_c)));
for i = 2:n_LingValsIn-1
    u_x1(i) = exp(-((x1-c_IN(i)).^2)/(2*sigma_IN^2))';
end
u_x1(n_LingValsIn) = 1./(1+exp(-sigmoid_a*(x1-sigmoid_c)));

u_x2 = zeros(1,n_LingValsIn);
u_x2(1) = 1./(1+exp(sigmoid_a*(x2+sigmoid_c)));
for i = 2:n_LingValsIn-1
    u_x2(i) = exp(-((x2-c_IN(i)).^2)/(2*sigma_IN^2))';
end
u_x2(n_LingValsIn) = 1./(1+exp(-sigmoid_a*(x2-sigmoid_c)));

% =========================================================================
%% Fuzzy Inference Mechanism 

if implied_flag==1  %Use Implied Fuzzy Sets (Eq. 4.49, 4.50)
    area = 0;
    num  = 0;
    den  = 0;
    uA   = 0;
    
    for i=1:n_LingValsIn
        for j=1:n_LingValsIn
            
            % Antecedent Composition (MIN)
            uA = min(u_x2(i),u_x1(j));
            
            % Implication (MIN Operator - Mamdani)
            implied_mf = min(out_mf(rules(i,j),:),uA);
            
            %Plotting Implied Membership Functions
            %{
            figure(1+i)
            clf
            plot(implied_mf)
            title('Implied Membership Functions')
            %}
            
            % Find area under all implied fuzzy sets
            area = sum( (1/(n_pts))*(implied_mf(1:n_pts)));
            
            % Defuzzification (COG Numerator)
            num = num + cMISO(rules(i,j))*area;
            
            % Defuzzification (COG Denominator)
            den = den + area;
        end
    end
    
    % Crisp Output
    lambda = num/den;
    
% =========================================================================    
else %Use Aggregated Implied Fuzzy Sets

    for i=1:n_LingValsIn
        for j=1:n_LingValsIn
            
            % Antecedent Composition (MIN)
            uA = min(u_x2(i),u_x1(j));
            
            % Implication (MIN - Mamdani)
            implied_mf((i-1)*9+j,1:n_pts) = min(out_mf(rules(i,j),:),uA);
        end
    end
    
    % Aggregation (MAX)
    for i=1:n_pts
        aggregate_mf(1,i)= max(implied_mf(:,i));
    end
    
    % Defuzzification (COA Numerator)
    Y_aggregate_mf = Y_pts.*aggregate_mf;
    num = sum((1/n_pts)*(Y_aggregate_mf(1:n_pts)));
    
    % Defuzzification (COA Denominator)
    den = sum((1/n_pts)*(aggregate_mf(1:n_pts)));
    
    % Crisp Output
    lambda = num/den;    
end

end
% =========================================================================