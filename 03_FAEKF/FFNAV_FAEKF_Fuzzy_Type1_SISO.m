function [lambda] = FFNAV_FAEKF_Fuzzy_Type1_SISO(x, FUZZY_props)
% FFNAV Fuzzy Controller Type 1 - SISO ====================================
% Description: This function completes the fuzzification, implication, and
% defuzzification step of the FLS. Given a single normalized input, the 
% function returns the corresponding normalized output of the fuzzy system.
%
% Inputs:
%   x            - Normalized input to the fuzzy system (one input) 
%   FUZZY_props. - Structure containing properties of the FLS
%
% Outputs:
%   lambda       - Normalized output of the fuzzy system
%
% Other Function Calls:
%   None
%
% Created by:  Cory Fraser - NOV 15, 2017
% Latest Edit: Cory Fraser - JUL 04, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% Initialize Parameters
rules_SISO      = FUZZY_props.rules_SISO;
n_pts           = FUZZY_props.n_pts;
X_pts           = FUZZY_props.X_pts;
Y_pts           = FUZZY_props.Y_pts;
out_mf_SISO     = FUZZY_props.out_mf_SISO;
implied_mf_SISO = FUZZY_props.implied_mf_SISO;
aggregate_mf    = FUZZY_props.aggregate_mf;
c_IN            = FUZZY_props.c_IN;
c_SISO          = FUZZY_props.c_SISO;
implied_flag    = FUZZY_props.implied_flag;
n_LingValsIn    = FUZZY_props.n_LingValsIn;
n_LingValsOut   = FUZZY_props.n_LingValsOut_SISO;
sigma_IN        = FUZZY_props.sigma_IN;
sigmoid_a       = FUZZY_props.sigmoid_a     ;
sigmoid_c       = FUZZY_props.sigmoid_c;

% =========================================================================
%% Fuzzification (Eq. 4.45, 4.46)

u_x1 = zeros(1,n_LingValsIn);

u_x1(1) = 1./(1+exp(sigmoid_a*(x+sigmoid_c)));
for i = 2:n_LingValsIn-1
    u_x1(i) = exp(-((x-c_IN(i)).^2)/(2*sigma_IN^2))';
end
u_x1(n_LingValsIn) = 1./(1+exp(-sigmoid_a*(x-sigmoid_c)));

% =========================================================================
%% Fuzzy Inference Mechanism 

if implied_flag == 1   %Use Implied Fuzzy Sets (Eq. 4.49, 4.50)
    area = 0;
    num = 0;
    den = 0;    
    for i=1:n_LingValsIn
        
        % Implication (MIN Operator - Mamdani)
        implied_mf_SISO = min(  out_mf_SISO(rules_SISO(i),:),  u_x1(i)  );
        
        %Plotting Implied Membership Functions
        %{
        figure(1+i)
        clf
        plot(implied_mfSISO)
        title('Implied Membership Functions')
        %}
        % Find area under all implied fuzzy sets
        area = sum( (1/(n_pts))*(implied_mf_SISO(1:n_pts)));
        
        % Defuzzification (COG Numerator)
        num = num + c_SISO(rules_SISO(i))*area;
        
        % Defuzzification (COG Denominator)
        den = den + area;
    end
    
    % Crisp Output
    lambda = num/den;
    
% =========================================================================
else % Use Aggregated Implied Fuzzy Set
    
    for i=1:n_LingValsOut
        
        % Implication (MIN - Mamdani)
        implied_mf_SISO(i,1:n_pts) = min(out_mf_SISO(rules_SISO(i),:),u_x1(i));
    
        %Plotting Output Membership Functions
        %{
        figure(1+i)
        %plot(implied_mfSISO(i,:))
        %title('Implied Membership Functions')
        %}
    end
        
    % Aggregation (MAX)
    for i=1:n_pts
        aggregate_mf(1,i) = max(implied_mf_SISO(:,i));
    end
    
    %Plotting Aggregate Membership Function
    %{
    figure(n_LingValsOut + 2)
    clf
    plot(aggregate_mf)
    title('Aggregate Membership Function')
    %}
    
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