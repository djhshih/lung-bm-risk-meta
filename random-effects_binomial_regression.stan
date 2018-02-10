data {
	int<lower=1> I;       // number of studies
	int<lower=1> J;       // number of covariates

	int<lower=0> m[I];    // numbers of successes
	int<lower=1> n[I];    // numbers of trials
	real x[I, J];         // covariates
}

parameters {
	real mu;              // overall effect
	real<lower=0> tau;    // between-study sd
	real beta[J];

	real z[I];
}

transformed parameters {
	real theta[I];          // study mean effect

	for (i in 1:I) {
		real t;
		t = mu + tau*z[i];
		for (j in 1:J) {
			t += beta[j] * x[i,j];
		}
		theta[i] = t;
	}
	
}

model {
	mu ~ normal(0, 5);
	tau ~ cauchy(0, 5);
	beta ~ normal(0, 5);
	z ~ normal(0, 1);
	m ~ binomial_logit(n, theta);
}

