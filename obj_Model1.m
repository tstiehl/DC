%% cost function for Model 1

function residuum = obj_Model1(par)

    global PC_SI_totalDC_24 PC_SI_totalDC_48 PC_SI_totalDC_72 % handle experimental data as global variable

    p = par(1); % proliferation rate
    mu = par(2); % outflux rate

    % calculate percent of photo converted cells at times 1 day, 2 days, 3 days
    converted_percent_model = photoconverted_Model1([1,2,3], p, mu);

    % calucalte residuals (weighted by variance of the data)
    r1 =  (PC_SI_totalDC_24 - converted_percent_model(1)) * (PC_SI_totalDC_24 - converted_percent_model(1))' / var(PC_SI_totalDC_24);
    r2 =  (PC_SI_totalDC_48 - converted_percent_model(2)) * (PC_SI_totalDC_48 - converted_percent_model(2))' / var(PC_SI_totalDC_48);
    r3 =  (PC_SI_totalDC_72 - converted_percent_model(3)) * (PC_SI_totalDC_72 - converted_percent_model(3))' / var(PC_SI_totalDC_72);

    residuum =  r1+r2+r3;




