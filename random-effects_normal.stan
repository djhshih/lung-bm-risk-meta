data {
	int<lower=1> I;       // number of studies
	real y[I];            // observed outcome
	real<lower=0> w[I];   // observation sd
}

parameters {
	real mu;              // overall effect
	real<lower=0> tau;    // between-study sd

	real z[I];
}

transformed parameters {
	real theta[I];        // study mean effect

	for (i in 1:I) {
		theta[i] = mu + tau * z[i];
	}
}

model {
	mu ~ normal(0, 10);
	tau ~ cauchy(0, 2);
	z ~ normal(0, 1);
	y ~ normal(theta, w);
}

