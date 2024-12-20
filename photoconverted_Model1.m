%% Calculates percent of photo-converted cell in small intestine according to Model 1 at time t
% mu: outflux rate
% p: proliferation rate

function percent = photoconverted_Model1(t, p, mu)

    percent = exp((p-mu)*t)*100; % equation (13) of Modeling Supplement 
