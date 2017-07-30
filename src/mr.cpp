#include "lib/calculate_hubs.h"
#include <Rcpp.h>

class r_matrix : public MatrixBase<double> {
public:
  explicit r_matrix(const Rcpp::NumericMatrix &data)
    : nrow(data.rows()), ncol(data.cols()) {
    ndata.resize(nrow);
    for (std::size_t index = 0; index < nrow; index++)
      ndata[index] = Rcpp::as<std::vector<double>>(
        static_cast<Rcpp::NumericVector>(data(index, Rcpp::_)));
  }
  
  std::size_t row_size() const override { return nrow; }
  std::size_t col_size() const override { return ncol; }
  
  const double *operator[](std::size_t index) const override {
    return &*ndata[index].begin();
  }
  
private:
  const std::size_t nrow;
  const std::size_t ncol;
  std::vector<std::vector<double>> ndata;
};

namespace Rcpp {
template <> SEXP wrap(const HubInfo &x) {
  return Rcpp::wrap(Rcpp::NumericVector::create(
      Rcpp::Named("if_hub") = x.if_hub,
      Rcpp::Named("cluster_bound") = x.cluster_bound,
      Rcpp::Named("index") = x.index));
};
} // namespace Rcpp

// [[Rcpp::export]]
Rcpp::List GetHubInfos(Rcpp::NumericMatrix m, double tn_p, std::size_t k) {
  std::vector<HubInfo> hubs;
  std::vector<std::size_t> mr_id;
  std::vector<std::vector<int>> mr_em;
  std::tie(hubs, mr_id, mr_em) = calculate_hubs(r_matrix(m), tn_p, k);
  Rcpp::List z = Rcpp::List::create(hubs, mr_id, mr_em);
  return z;
}

// [[Rcpp::export]]
Rcpp::List GetMrEm(Rcpp::NumericMatrix m, double tn_p, std::size_t k) {
  std::vector<std::size_t> mr_id;
  std::vector<std::vector<int>> mr_em;
  std::tie(mr_em, mr_id) = get_mr_em(r_matrix(m), tn_p, k);
  Rcpp::List z = Rcpp::List::create(mr_em, mr_id);
  return z;
}

// [[Rcpp::export]]
Rcpp::List GetMrhcaFast(Rcpp::NumericMatrix m, std::vector<std::vector<int>> mr_em, std::size_t k, std::size_t step_size = 50) {
  std::vector<HubInfo> hubs = MRHCA_fast_version(r_matrix(m), mr_em, k, step_size);
  return Rcpp::List::create(hubs);
}