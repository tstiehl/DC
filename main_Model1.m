%% Fit Model 1 to photo-conversion data

clear all
close all
global PC_SI_totalDC_24 PC_SI_totalDC_48 PC_SI_totalDC_72

% load photo-conversion data
dataPC;

% boundaries for parameter estimation
% first parameter: p
% second parameter: mu
lower_bd_LH = [0.1, 0];
upper_bd_LH = [0.1, 10];
% model only depends on mu-p, therefore, fix p



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

fun = @(x)obj_Model1(x); % set cost function 
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
% save('fit_Model1.mat', 'LH', 'initguessmatrix', 'params','fit_res', 'parmatrix', 'best_fit', "grad_all", "lambda_all","hessian_all")


p = best_fit(1); % fitted proliferation rate
mu = best_fit(2); % fitted outflux rate

% plot best fit
T_pc = linspace(0,15,200); % time points for plotting
converted_percent_model = photoconverted_Model1(T_pc, p, mu); % calculated percent of photo-converted cells

figure()
hold on
plot(T_pc, converted_percent_model,'r-','LineWidth',3) % plot fitted model
plot(1,PC_SI_totalDC_24,'bo','LineWidth',3) % plot data t = 1 day
plot(2,PC_SI_totalDC_48,'bo','LineWidth',3) % plot data t = 2 days
plot(3,PC_SI_totalDC_72,'bo','LineWidth',3) % plot data t = 3 days
influx = mu - p; % fraction of population entering small intestine per day

title(['p: ', num2str(p), '  influx: ', num2str(influx),' outflux: ',num2str(mu) ,'  half-life: ', num2str(log(2)/(mu-p)), '  divisions: ',num2str(log(mu/influx)/log(2))])
%                                                                                                                                           equation (3)
% rates have unit 1/day, half-life has unit days
ylabel('% photoconverted DCs')
xlabel('time [days]')
% saveas(gcf,['Model1','.png'])
% saveas(gcf,['Model1','.svg'])

