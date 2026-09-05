// Copyright (C) 2026 JLinAlg contributors
// SPDX-License-Identifier: GPL-2.0-or-later
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
    DATA_VECTOR(y);
    DATA_MATRIX(X_count);
    DATA_MATRIX(X_zero);
    DATA_SPARSE_MATRIX(Z);
    DATA_SPARSE_MATRIX(A_inverse);
    DATA_SCALAR(log_determinant_A_inverse);
    PARAMETER_VECTOR(beta_count);
    PARAMETER_VECTOR(beta_zero);
    PARAMETER_VECTOR(a_count);
    PARAMETER_VECTOR(a_zero);
    PARAMETER(log_sd_count);
    PARAMETER(log_sd_zero);
    PARAMETER(fisher_z);

    Type rho = Type(0.99) * tanh(fisher_z);
    Type denominator = Type(1.0) - rho * rho;
    Type variance_count = exp(Type(2.0) * log_sd_count);
    Type variance_zero = exp(Type(2.0) * log_sd_zero);
    vector<Type> qa_count = A_inverse * a_count;
    vector<Type> qa_zero = A_inverse * a_zero;
    Type quadratic = (a_count * qa_count).sum() / variance_count
        + (a_zero * qa_zero).sum() / variance_zero
        - Type(2.0) * rho * (a_count * qa_zero).sum()
            / exp(log_sd_count + log_sd_zero);
    Type nll = Type(0.5) * quadratic / denominator
        - log_determinant_A_inverse
        + Type(a_count.size()) * (log_sd_count + log_sd_zero)
        + Type(0.5 * a_count.size()) * log(denominator)
        + Type(a_count.size()) * log(Type(2.0 * M_PI));

    vector<Type> eta_count = X_count * beta_count + Z * a_count;
    vector<Type> eta_zero = X_zero * beta_zero + Z * a_zero;
    for (int row = 0; row < y.size(); row++) {
        Type mean = exp(eta_count(row));
        Type pi = invlogit(eta_zero(row));
        if (y(row) == Type(0.0))
            nll -= log(pi + (Type(1.0) - pi) * exp(-mean));
        else
            nll -= log(Type(1.0) - pi) + dpois(y(row), mean, true);
    }
    ADREPORT(beta_count);
    ADREPORT(beta_zero);
    ADREPORT(log_sd_count);
    ADREPORT(log_sd_zero);
    ADREPORT(rho);
    return nll;
}
