%% Calculates percent of photo-converted cell in small intestine according to Model 3 at time t
% p: parameter quantifying proliferation rate, see equation (57) of the Modeling Supplement
% T: time DCs spend in small intestine (denoted by $\hat \tau$ in the Modeling Supplement)
% gamma: premature exit rate

function percent = photoconverted_Model3(t, p, T, gamma)
    t = min(t, T);
    fun = @(tau) exp(tau.^2/2*p-gamma*tau);
    percent = (1 - integral(fun,0,t)/integral(fun,0,T))*100; % equation (77) of Modeling Supplement
