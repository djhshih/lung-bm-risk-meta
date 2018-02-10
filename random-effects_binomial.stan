data {
	int<lower=1> I;       // number of studies
	int<lower=0> m[I];    // number of successes
	int<lower=1> n[I];    // number of trials
}

parameters {
	real mu;              // overall effect
	real<lower=0> tau;    // between-study sd

	real z[I];
}

transformed parameters {
	real theta[I];          // study mean effect

	for (i in 1:I) {
		theta[i] = mu + tau * z[i];
	}
	
}

model {
	mu ~ normal(0, 10);
	tau ~ cauchy(0, 5);
	z ~ normal(0, 1);
	m ~ binomial_logit(n, theta);
}

