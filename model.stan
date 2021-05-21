functions{
    // discretised truncated lognormal pmf
    vector discretised_lognormal_pmf(int[] y, real mu, real sigma, int max_val) {
    int n = num_elements(y);
    vector[n] pmf;
    real small = 1e-5;
    real c_sigma = sigma < small ? small : sigma;
    real c_mu = mu < small ? small : mu;
    vector[n] adj_y = to_vector(y) + small;
    vector[n] upper_y = (log(adj_y + 1) - c_mu) / c_sigma;
    vector[n] lower_y = (log(adj_y) - c_mu) / c_sigma;
    real max_cdf = normal_cdf((log(max_val + small) - c_mu) / c_sigma, 0.0, 1.0);
    real min_cdf = normal_cdf((log(small) - c_mu) / c_sigma, 0.0, 1.0);
    real trunc_cdf = max_cdf - min_cdf;
    for (i in 1:n) {
        pmf[i] = (normal_cdf(upper_y[i], 0.0, 1.0) - normal_cdf(lower_y[i], 0.0, 1.0)) /
        trunc_cdf;
    }
    return(pmf);
    }

}
data {
    int t;
    int cases_obs[t];
    int cases_seq[t];
    int cases_b1672[t];
    int imp_b1672[t];
}

parameters {
    real <lower = 0> R;
    real <lower = 0> b1672_mod;
    real <lower = 0> imp_mod;
    real si_logmean;
    real <lower = 0> si_logsd;
    real <lower = 0> recip_phi;
}

transfomed parameters {
    vector[t] imp_linked_b1672;
    vector[t] exp_b1672;
    vector[t] exp_b117;
    vector si[30];
    real phi;

    // discretised serial interval
    si = discretised_lognormal_pmf(1:30, si_logmean, si_logsd, 30);

    // b1672 cases directly driven by imports
    for (s in 1:t) {    
        int index = min(t - s, si_t);
        imp_linked_b1672[s:index] = imp_linked_b1672[s:index] + 
            imp_b1672[s] * si[1:index] * R * b1672_mod + imp_mod;
    }

    // b1672 cases from transmission
    for (s in 1:t) {    
        int index = min(t - s, si_t);
        exp_b1672[s:index] = exp_b1672[s:index] + 
            exp_b1672[s] * si[1:index] * R * b1672_mod;
    }

    // b117 cases from transmission
    exp_b117 = cases_obs[1];
    for (s in 2:t) {
        int index = min(t - s, si_t);
        exp_b117[s:index] = exp_b117[s:index] + 
            exp_b117[s] * si[1:index] * R;
    }
    // convert overdispersion to correct scale
    phi = 1 ./ sqrt(recip_phi);

}

model {
    si_logmean ~ normal(1.65, 0.1);
    si_logsd ~ normal(0.273, 0.05) T[0, ];

    R ~ lognormal(0, 0.25);
    imp_mod ~ lognormal(0, 1); 
    b1672_mod ~ lognormal(0, 1);

    cases_b1672 ~ hypergeometric(cases_seq, exp_b1672, cases_obs - exp_b1672);
    cases_obs ~ neg_binomial_2(exp_b1672 + exp_b117, phi);
}

generated quantities {
   
}