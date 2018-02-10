// Mixed meta-regression with errors-in-variable model extension

data {
	int<lower=1> I;       // number of studies
	int<lower=1> K1;       // number of groups
	int<lower=1> K2;       // number of groups

	int<lower=0> m[I];    // numbers of successes
	int<lower=1> n[I];    // numbers of trials

	vector<lower=0>[K1] alpha1;      // prior counts for each group
	vector<lower=0>[K2] alpha2;      // prior counts for each group
	int<lower=0> X1[I, K1];         // group counts
	int<lower=0> X2[I, K2];         // group counts
}

parameters {
	real mu;              // overall effect
	real<lower=0> tau;    // between-study sd
	real beta1[K1-1];
	real beta2[K2-1];

	simplex[K1] pi1[I];           // covariate level probabilities
	simplex[K2] pi2[I];           // covariate level probabilities

	real z[I];
}

transformed parameters {
	real theta[I];          // study mean effect

	for (i in 1:I) {
		real t;
		t = mu + tau*z[i];
		for (k in 1:(K1-1)) {
			t += beta1[k] * logit(pi1[i, k+1]);
		}
		for (k in 1:(K2-1)) {
			t += beta2[k] * logit(pi2[i, k+1]);
		}
		theta[i] = t;
	}
	
}

model {
	mu ~ normal(0, 5);
	tau ~ cauchy(0, 5);
	beta1 ~ normal(0, 5);
	beta2 ~ normal(0, 5);
	for (i in 1:I) {
		pi1[i] ~ dirichlet(alpha1);
		pi2[i] ~ dirichlet(alpha2);
		X1[i] ~ multinomial(pi1[i]);
		X2[i] ~ multinomial(pi2[i]);
	}
	z ~ normal(0, 1);
	m ~ binomial_logit(n, theta);
}

