// Generated by rstantools.  Do not edit by hand.

/*
    pkg is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    pkg is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with pkg.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.19.1
#include <stan/model/model_header.hpp>
namespace model_bayes_cpm_mod_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_bayes_cpm_mod");
    reader.add_event(7, 7, "include", "/sub/bayes_cpm_funs.stan");
    reader.add_event(7, 0, "start", "/sub/bayes_cpm_funs.stan");
    reader.add_event(86, 79, "end", "/sub/bayes_cpm_funs.stan");
    reader.add_event(86, 8, "restart", "model_bayes_cpm_mod");
    reader.add_event(131, 51, "end", "model_bayes_cpm_mod");
    return reader;
}
template <typename T0__>
typename boost::math::tools::promote_args<T0__>::type
CDF_polr(const T0__& x,
             const int& link, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        current_statement_begin__ = 12;
        if (as_bool(logical_eq(link, 1))) {
            current_statement_begin__ = 12;
            return stan::math::promote_scalar<fun_return_scalar_t__>(inv_logit(x));
        } else if (as_bool(logical_eq(link, 2))) {
            current_statement_begin__ = 13;
            return stan::math::promote_scalar<fun_return_scalar_t__>(Phi(x));
        } else if (as_bool(logical_eq(link, 3))) {
            current_statement_begin__ = 14;
            return stan::math::promote_scalar<fun_return_scalar_t__>(gumbel_cdf(x, 0, 1));
        } else if (as_bool(logical_eq(link, 4))) {
            current_statement_begin__ = 15;
            return stan::math::promote_scalar<fun_return_scalar_t__>(inv_cloglog(x));
        } else if (as_bool(logical_eq(link, 5))) {
            current_statement_begin__ = 16;
            return stan::math::promote_scalar<fun_return_scalar_t__>(cauchy_cdf(x, 0, 1));
        } else {
            current_statement_begin__ = 17;
            std::stringstream errmsg_stream__;
            errmsg_stream__ << "Invalid link";
            throw std::domain_error(errmsg_stream__.str());
        }
        current_statement_begin__ = 18;
        return stan::math::promote_scalar<fun_return_scalar_t__>(x);
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct CDF_polr_functor__ {
    template <typename T0__>
        typename boost::math::tools::promote_args<T0__>::type
    operator()(const T0__& x,
             const int& link, std::ostream* pstream__) const {
        return CDF_polr(x, link, pstream__);
    }
};
template <typename T0__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic, 1>
make_cutpoints(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& probabilities,
                   const int& ncat,
                   const int& link, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 29;
        validate_non_negative_index("cutpoints", "(ncat - 1)", (ncat - 1));
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> cutpoints((ncat - 1));
        stan::math::initialize(cutpoints, DUMMY_VAR__);
        stan::math::fill(cutpoints, DUMMY_VAR__);
        current_statement_begin__ = 30;
        local_scalar_t__ running_sum(DUMMY_VAR__);
        (void) running_sum;  // dummy to suppress unused var warning
        stan::math::initialize(running_sum, DUMMY_VAR__);
        stan::math::fill(running_sum, DUMMY_VAR__);
        stan::math::assign(running_sum,0);
        current_statement_begin__ = 31;
        if (as_bool(logical_eq(link, 1))) {
            current_statement_begin__ = 31;
            for (int i = 1; i <= (ncat - 1); ++i) {
                current_statement_begin__ = 32;
                stan::math::assign(running_sum, (running_sum + get_base1(probabilities, i, "probabilities", 1)));
                current_statement_begin__ = 33;
                stan::model::assign(cutpoints, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            logit(running_sum), 
                            "assigning variable cutpoints");
            }
        } else if (as_bool(logical_eq(link, 2))) {
            current_statement_begin__ = 35;
            for (int i = 1; i <= (ncat - 1); ++i) {
                current_statement_begin__ = 36;
                stan::math::assign(running_sum, (running_sum + get_base1(probabilities, i, "probabilities", 1)));
                current_statement_begin__ = 37;
                stan::model::assign(cutpoints, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            inv_Phi(running_sum), 
                            "assigning variable cutpoints");
            }
        } else if (as_bool(logical_eq(link, 3))) {
            current_statement_begin__ = 39;
            for (int i = 1; i <= (ncat - 1); ++i) {
                current_statement_begin__ = 40;
                stan::math::assign(running_sum, (running_sum + get_base1(probabilities, i, "probabilities", 1)));
                current_statement_begin__ = 41;
                stan::model::assign(cutpoints, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            -(stan::math::log(-(stan::math::log(running_sum)))), 
                            "assigning variable cutpoints");
            }
        } else if (as_bool(logical_eq(link, 4))) {
            current_statement_begin__ = 43;
            for (int i = 1; i <= (ncat - 1); ++i) {
                current_statement_begin__ = 44;
                stan::math::assign(running_sum, (running_sum + get_base1(probabilities, i, "probabilities", 1)));
                current_statement_begin__ = 45;
                stan::model::assign(cutpoints, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            stan::math::log(-(log1m(running_sum))), 
                            "assigning variable cutpoints");
            }
        } else if (as_bool(logical_eq(link, 5))) {
            current_statement_begin__ = 47;
            for (int i = 1; i <= (ncat - 1); ++i) {
                current_statement_begin__ = 48;
                stan::math::assign(running_sum, (running_sum + get_base1(probabilities, i, "probabilities", 1)));
                current_statement_begin__ = 49;
                stan::model::assign(cutpoints, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            stan::math::tan((stan::math::pi() * (running_sum - 0.5))), 
                            "assigning variable cutpoints");
            }
        } else {
            current_statement_begin__ = 51;
            std::stringstream errmsg_stream__;
            errmsg_stream__ << "invalid link";
            throw std::domain_error(errmsg_stream__.str());
        }
        current_statement_begin__ = 52;
        return stan::math::promote_scalar<fun_return_scalar_t__>(cutpoints);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct make_cutpoints_functor__ {
    template <typename T0__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic, 1>
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& probabilities,
                   const int& ncat,
                   const int& link, std::ostream* pstream__) const {
        return make_cutpoints(probabilities, ncat, link, pstream__);
    }
};
template <typename T2__, typename T4__>
Eigen::Matrix<typename boost::math::tools::promote_args<T2__, T4__>::type, Eigen::Dynamic, 1>
loglik(const std::vector<int>& Ylev,
           const int& N,
           const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& cutpoints,
           const int& ncat,
           const Eigen::Matrix<T4__, Eigen::Dynamic, 1>& eta,
           const int& link, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T2__, T4__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 58;
        validate_non_negative_index("ll", "N", N);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> ll(N);
        stan::math::initialize(ll, DUMMY_VAR__);
        stan::math::fill(ll, DUMMY_VAR__);
        current_statement_begin__ = 59;
        for (int n = 1; n <= N; ++n) {
            current_statement_begin__ = 60;
            if (as_bool(logical_eq(get_base1(Ylev, n, "Ylev", 1), 1))) {
                current_statement_begin__ = 61;
                stan::model::assign(ll, 
                            stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                            stan::math::log(CDF_polr((get_base1(cutpoints, 1, "cutpoints", 1) - get_base1(eta, n, "eta", 1)), link, pstream__)), 
                            "assigning variable ll");
            } else if (as_bool(logical_eq(get_base1(Ylev, n, "Ylev", 1), ncat))) {
                current_statement_begin__ = 64;
                stan::model::assign(ll, 
                            stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                            log1m(CDF_polr((get_base1(cutpoints, (ncat - 1), "cutpoints", 1) - get_base1(eta, n, "eta", 1)), link, pstream__)), 
                            "assigning variable ll");
            } else {
                current_statement_begin__ = 67;
                stan::model::assign(ll, 
                            stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                            stan::math::log((CDF_polr((get_base1(cutpoints, get_base1(Ylev, n, "Ylev", 1), "cutpoints", 1) - get_base1(eta, n, "eta", 1)), link, pstream__) - CDF_polr((get_base1(cutpoints, (get_base1(Ylev, n, "Ylev", 1) - 1), "cutpoints", 1) - get_base1(eta, n, "eta", 1)), link, pstream__))), 
                            "assigning variable ll");
            }
        }
        current_statement_begin__ = 69;
        return stan::math::promote_scalar<fun_return_scalar_t__>(ll);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct loglik_functor__ {
    template <typename T2__, typename T4__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T2__, T4__>::type, Eigen::Dynamic, 1>
    operator()(const std::vector<int>& Ylev,
           const int& N,
           const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& cutpoints,
           const int& ncat,
           const Eigen::Matrix<T4__, Eigen::Dynamic, 1>& eta,
           const int& link, std::ostream* pstream__) const {
        return loglik(Ylev, N, cutpoints, ncat, eta, link, pstream__);
    }
};
template <typename T2__, typename T4__, typename T5__>
Eigen::Matrix<typename boost::math::tools::promote_args<T2__, T4__, T5__>::type, Eigen::Dynamic, 1>
loglik_hier(const std::vector<int>& Ylev,
                const int& N,
                const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& cutpoints,
                const int& ncat,
                const Eigen::Matrix<T4__, Eigen::Dynamic, 1>& eta,
                const Eigen::Matrix<T5__, Eigen::Dynamic, 1>& zeta,
                const std::vector<int>& jj,
                const int& link, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T2__, T4__, T5__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 74;
        validate_non_negative_index("ll", "N", N);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> ll(N);
        stan::math::initialize(ll, DUMMY_VAR__);
        stan::math::fill(ll, DUMMY_VAR__);
        current_statement_begin__ = 75;
        for (int n = 1; n <= N; ++n) {
            current_statement_begin__ = 76;
            if (as_bool(logical_eq(get_base1(Ylev, n, "Ylev", 1), 1))) {
                current_statement_begin__ = 77;
                stan::model::assign(ll, 
                            stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                            stan::math::log(CDF_polr(((get_base1(cutpoints, 1, "cutpoints", 1) - get_base1(eta, n, "eta", 1)) - get_base1(zeta, get_base1(jj, n, "jj", 1), "zeta", 1)), link, pstream__)), 
                            "assigning variable ll");
            } else if (as_bool(logical_eq(get_base1(Ylev, n, "Ylev", 1), ncat))) {
                current_statement_begin__ = 80;
                stan::model::assign(ll, 
                            stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                            stan::math::log((1 - CDF_polr(((get_base1(cutpoints, (ncat - 1), "cutpoints", 1) - get_base1(eta, n, "eta", 1)) - get_base1(zeta, get_base1(jj, n, "jj", 1), "zeta", 1)), link, pstream__))), 
                            "assigning variable ll");
            } else {
                current_statement_begin__ = 83;
                stan::model::assign(ll, 
                            stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                            stan::math::log((CDF_polr(((get_base1(cutpoints, get_base1(Ylev, n, "Ylev", 1), "cutpoints", 1) - get_base1(eta, n, "eta", 1)) - get_base1(zeta, get_base1(jj, n, "jj", 1), "zeta", 1)), link, pstream__) - CDF_polr(((get_base1(cutpoints, (get_base1(Ylev, n, "Ylev", 1) - 1), "cutpoints", 1) - get_base1(eta, n, "eta", 1)) - get_base1(zeta, get_base1(jj, n, "jj", 1), "zeta", 1)), link, pstream__))), 
                            "assigning variable ll");
            }
        }
        current_statement_begin__ = 85;
        return stan::math::promote_scalar<fun_return_scalar_t__>(ll);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct loglik_hier_functor__ {
    template <typename T2__, typename T4__, typename T5__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T2__, T4__, T5__>::type, Eigen::Dynamic, 1>
    operator()(const std::vector<int>& Ylev,
                const int& N,
                const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& cutpoints,
                const int& ncat,
                const Eigen::Matrix<T4__, Eigen::Dynamic, 1>& eta,
                const Eigen::Matrix<T5__, Eigen::Dynamic, 1>& zeta,
                const std::vector<int>& jj,
                const int& link, std::ostream* pstream__) const {
        return loglik_hier(Ylev, N, cutpoints, ncat, eta, zeta, jj, link, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_bayes_cpm_mod : public prob_grad {
private:
        int N;
        int ncat;
        std::vector<int> Ylev;
        int link;
        int K;
        matrix_d Q;
        double alpha;
public:
    model_bayes_cpm_mod(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_bayes_cpm_mod(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_bayes_cpm_mod_namespace::model_bayes_cpm_mod";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 89;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 90;
            context__.validate_dims("data initialization", "ncat", "int", context__.to_vec());
            ncat = int(0);
            vals_i__ = context__.vals_i("ncat");
            pos__ = 0;
            ncat = vals_i__[pos__++];
            current_statement_begin__ = 91;
            validate_non_negative_index("Ylev", "N", N);
            context__.validate_dims("data initialization", "Ylev", "int", context__.to_vec(N));
            Ylev = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("Ylev");
            pos__ = 0;
            size_t Ylev_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < Ylev_k_0_max__; ++k_0__) {
                Ylev[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 92;
            context__.validate_dims("data initialization", "link", "int", context__.to_vec());
            link = int(0);
            vals_i__ = context__.vals_i("link");
            pos__ = 0;
            link = vals_i__[pos__++];
            current_statement_begin__ = 93;
            context__.validate_dims("data initialization", "K", "int", context__.to_vec());
            K = int(0);
            vals_i__ = context__.vals_i("K");
            pos__ = 0;
            K = vals_i__[pos__++];
            current_statement_begin__ = 94;
            validate_non_negative_index("Q", "N", N);
            validate_non_negative_index("Q", "K", K);
            context__.validate_dims("data initialization", "Q", "matrix_d", context__.to_vec(N,K));
            Q = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N, K);
            vals_r__ = context__.vals_r("Q");
            pos__ = 0;
            size_t Q_j_2_max__ = K;
            size_t Q_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < Q_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < Q_j_1_max__; ++j_1__) {
                    Q(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 95;
            context__.validate_dims("data initialization", "alpha", "double", context__.to_vec());
            alpha = double(0);
            vals_r__ = context__.vals_r("alpha");
            pos__ = 0;
            alpha = vals_r__[pos__++];
            check_greater_or_equal(function__, "alpha", alpha, 0);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 99;
            validate_non_negative_index("pi", "ncat", ncat);
            num_params_r__ += (ncat - 1);
            current_statement_begin__ = 100;
            validate_non_negative_index("b", "K", K);
            num_params_r__ += K;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_bayes_cpm_mod() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 99;
        if (!(context__.contains_r("pi")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable pi missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("pi");
        pos__ = 0U;
        validate_non_negative_index("pi", "ncat", ncat);
        context__.validate_dims("parameter initialization", "pi", "vector_d", context__.to_vec(ncat));
        Eigen::Matrix<double, Eigen::Dynamic, 1> pi(ncat);
        size_t pi_j_1_max__ = ncat;
        for (size_t j_1__ = 0; j_1__ < pi_j_1_max__; ++j_1__) {
            pi(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.simplex_unconstrain(pi);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable pi: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 100;
        if (!(context__.contains_r("b")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable b missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("b");
        pos__ = 0U;
        validate_non_negative_index("b", "K", K);
        context__.validate_dims("parameter initialization", "b", "vector_d", context__.to_vec(K));
        Eigen::Matrix<double, Eigen::Dynamic, 1> b(K);
        size_t b_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < b_j_1_max__; ++j_1__) {
            b(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(b);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable b: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 99;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> pi;
            (void) pi;  // dummy to suppress unused var warning
            if (jacobian__)
                pi = in__.simplex_constrain(ncat, lp__);
            else
                pi = in__.simplex_constrain(ncat);
            current_statement_begin__ = 100;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> b;
            (void) b;  // dummy to suppress unused var warning
            if (jacobian__)
                b = in__.vector_constrain(K, lp__);
            else
                b = in__.vector_constrain(K);
            // transformed parameters
            current_statement_begin__ = 104;
            validate_non_negative_index("cutpoints", "(ncat - 1)", (ncat - 1));
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> cutpoints((ncat - 1));
            stan::math::initialize(cutpoints, DUMMY_VAR__);
            stan::math::fill(cutpoints, DUMMY_VAR__);
            current_statement_begin__ = 105;
            validate_non_negative_index("log_lik", "N", N);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> log_lik(N);
            stan::math::initialize(log_lik, DUMMY_VAR__);
            stan::math::fill(log_lik, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 107;
            stan::math::assign(cutpoints, make_cutpoints(pi, ncat, link, pstream__));
            current_statement_begin__ = 108;
            stan::math::assign(log_lik, loglik(Ylev, N, cutpoints, ncat, multiply(Q, b), link, pstream__));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 104;
            size_t cutpoints_j_1_max__ = (ncat - 1);
            for (size_t j_1__ = 0; j_1__ < cutpoints_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(cutpoints(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: cutpoints" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable cutpoints: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 105;
            size_t log_lik_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(log_lik(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: log_lik" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable log_lik: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            current_statement_begin__ = 115;
            lp_accum__.add(dirichlet_log(pi, rep_vector(alpha, ncat)));
            current_statement_begin__ = 124;
            lp_accum__.add(log_lik);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("pi");
        names__.push_back("b");
        names__.push_back("cutpoints");
        names__.push_back("log_lik");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(ncat);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back((ncat - 1));
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_bayes_cpm_mod_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> pi = in__.simplex_constrain(ncat);
        size_t pi_j_1_max__ = ncat;
        for (size_t j_1__ = 0; j_1__ < pi_j_1_max__; ++j_1__) {
            vars__.push_back(pi(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> b = in__.vector_constrain(K);
        size_t b_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < b_j_1_max__; ++j_1__) {
            vars__.push_back(b(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 104;
            validate_non_negative_index("cutpoints", "(ncat - 1)", (ncat - 1));
            Eigen::Matrix<double, Eigen::Dynamic, 1> cutpoints((ncat - 1));
            stan::math::initialize(cutpoints, DUMMY_VAR__);
            stan::math::fill(cutpoints, DUMMY_VAR__);
            current_statement_begin__ = 105;
            validate_non_negative_index("log_lik", "N", N);
            Eigen::Matrix<double, Eigen::Dynamic, 1> log_lik(N);
            stan::math::initialize(log_lik, DUMMY_VAR__);
            stan::math::fill(log_lik, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 107;
            stan::math::assign(cutpoints, make_cutpoints(pi, ncat, link, pstream__));
            current_statement_begin__ = 108;
            stan::math::assign(log_lik, loglik(Ylev, N, cutpoints, ncat, multiply(Q, b), link, pstream__));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t cutpoints_j_1_max__ = (ncat - 1);
                for (size_t j_1__ = 0; j_1__ < cutpoints_j_1_max__; ++j_1__) {
                    vars__.push_back(cutpoints(j_1__));
                }
                size_t log_lik_j_1_max__ = N;
                for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
                    vars__.push_back(log_lik(j_1__));
                }
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    static std::string model_name() {
        return "model_bayes_cpm_mod";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t pi_j_1_max__ = ncat;
        for (size_t j_1__ = 0; j_1__ < pi_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "pi" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t b_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < b_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "b" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t cutpoints_j_1_max__ = (ncat - 1);
            for (size_t j_1__ = 0; j_1__ < cutpoints_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "cutpoints" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t log_lik_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "log_lik" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t pi_j_1_max__ = (ncat - 1);
        for (size_t j_1__ = 0; j_1__ < pi_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "pi" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t b_j_1_max__ = K;
        for (size_t j_1__ = 0; j_1__ < b_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "b" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t cutpoints_j_1_max__ = (ncat - 1);
            for (size_t j_1__ = 0; j_1__ < cutpoints_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "cutpoints" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t log_lik_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "log_lik" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_bayes_cpm_mod_namespace::model_bayes_cpm_mod stan_model;
namespace torsten {
namespace dsolve {
    template<typename... Args>
    inline auto pmx_ode_group_mpi_functor::operator()(Args&&... args) const {
        dummy_functor f; return f(std::forward<Args>(args)...);
    }
}
}
#endif
