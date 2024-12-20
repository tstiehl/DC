%% Fit Model 3 to photo-conversion data
%% Impose different values for p

clear all
close all
global PC_SI_totalDC_24 PC_SI_totalDC_48 PC_SI_totalDC_72

% load photo-conversion data
dataPC;

% Imposed values of the parameter p. For each value in the list, a model fit is performed
p_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3];

for i = 1:length(p_list)
    % boundaries for parameter estimation
    % first parameter: p
    % second parameter: T (in the supplement denoted by $\hat \tau$), i.e., time spent in small intestine
    % third parameter: gamma (premature exit rate)
    lower_bd_LH = [p_list(i), 0, 0];
    upper_bd_LH = [p_list(i), 10, 10];

    disp(['p = ', num2str(p_list(i))])

    % number of multistarts (set to 15 for demonstration purposes)
    n_multistart = 15; % 200;

    fit_res = ones(1,n_multistart)*Inf; % array to store fit results
    start = {}; % cellarray to store initial guesses
    params = {}; % cell array to store estimated parameters
    lambda_all = {}; % cell array to store fit results
    hessian_all = {}; % cell array to store fit results
    grad_all = {}; % cell array to store fit results
    parmatrix = []; % array to store fitted parameters for further analysis
    initguessmatrix =  [];  % array to store inital guesses for further analysis

    LH = lhsdesign(n_multistart,length(lower_bd_LH)); % sample initial guesses using lating hypercube
    LH_sz = size(LH);

    fun = @(x)obj_Model3(x); % set cost function
    opts = optimoptions('fmincon', 'MaxFunctionEvaluations',50000); % set options for fmincon


    for k = 1:LH_sz(1) % loop over multistarts
        disp(k)
        LH_sample = LH(k,:); % extract hypercube sample for inital guess
        initguess = lower_bd_LH + LH_sample.*(upper_bd_LH-lower_bd_LH); % calculate inital guess

        [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,initguess,[],[],[],[],lower_bd_LH,upper_bd_LH,[],opts); % run optimization

        % store fit result
        if exitflag>0
            [val,pos] = max(fit_res);
            parmatrix = [parmatrix;[x,fval,exitflag]];
            initguessmatrix = [initguessmatrix; initguess];
            grad_all{k} = grad;
            lambda_all{k} = lambda;
            hessian_all{k} = hessian;
            if fval<val
                fit_res(pos) = fval;
                params{pos} = x;
            end
        end
    end

    % extract best fit
    [val,pos] = min(fit_res);
    best_fit = params{pos};

    % % save results
    % save(['fit_Model3_loop_',num2str(i),'.mat'], 'LH', 'initguessmatrix', 'params','fit_res', 'parmatrix', 'best_fit', "grad_all", "lambda_all","hessian_all")


    p = best_fit(1); % fitted proliferation parameter
    T = best_fit(2); % fitted residence time
    gamma = best_fit(3); % fitted premature exit rate

    % plot best fit
    T_pc = linspace(0,15,200); % time points for plotting
    converted_percent_model = arrayfun(@(t) photoconverted_Model3(t, p, T,gamma),T_pc); % calculated percent of photo-converted cells

    figure()
    hold on
    plot(T_pc, converted_percent_model,'r-','LineWidth',3) % plot fitted model
    plot(1,PC_SI_totalDC_24,'bo','LineWidth',3) % plot data t = 1 day
    plot(2,PC_SI_totalDC_48,'bo','LineWidth',3) % plot data t = 2 days
    plot(3,PC_SI_totalDC_72,'bo','LineWidth',3) % plot data t = 3 days

    fun = @(tau) exp(tau.^2/2*p - gamma * tau);
    influx = 1 / integral(fun,0,T); % fraction of population entering small intestine per day, see p. 13 of Modeling Supplement
    mature_outflux = influx * fun(T); % fraction of population entering small intestine per day, see p. 12 of Modeling Supplement
    premature_outflux = influx * gamma * integral(fun,0,T);
    title({string(['p: ', num2str(p), '  influx: ', num2str(influx),' mat. out: ',num2str(mature_outflux), ' premat. out: ',num2str(premature_outflux)]) ,string(['stay: ', num2str(T), ' max.  div: ',num2str(p*T^2/log(2)/2)])}, 'FontSize',10)
    %                                                                                                                                      equation (78)
    % rates have unit 1/day, residence time has unit days
    ylabel('% photoconverted DCs')
    xlabel('time [days]')
    % saveas(gcf,['Model3_loop_',num2str(i),'.png'])
    % saveas(gcf,['Model3_loop_',num2str(i),'.svg'])
end
