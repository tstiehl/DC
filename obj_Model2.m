%% cost function for Model 2

function residuum = obj_Model2(par)

    global PC_SI_totalDC_24 PC_SI_totalDC_48 PC_SI_totalDC_72 % handle experimental data as global variable

    p = par(1); % parameter quantifying proliferation rate, see equation (57) of the Modeling Supplement
    T = par(2); % time DCs spend in small intestine (denoted by $\hat \tau$ in the Modeling Supplement)

    % calculate percent of photo converted cells at times 1 day, 2 days, 3 days

    converted_percent_model = arrayfun(@(t) photoconverted_Model2(t, p, T),[1,2,3]);

    % calucalte residuals (weighted by variance of the data)
    r1 =  (PC_SI_totalDC_24 - converted_percent_model(1)) * (PC_SI_totalDC_24 - converted_percent_model(1))' / var(PC_SI_totalDC_24);
    r2 =  (PC_SI_totalDC_48 - converted_percent_model(2)) * (PC_SI_totalDC_48 - converted_percent_model(2))' / var(PC_SI_totalDC_48);
    r3 =  (PC_SI_totalDC_72 - converted_percent_model(3)) * (PC_SI_totalDC_72 - converted_percent_model(3))' / var(PC_SI_totalDC_72);

    residuum =  r1+r2+r3;




