#include <Rcpp.h>
using namespace Rcpp;

//' Revert a treated-to-control distance list.
//'
//' @param n_t Number of treated units
//' @param n_c Number of control units
//' @param startn Vector of starting nodes of edges
//' @param endn Vector of ending nodes of edges
//' @param d Vector of cost associated with edges
//' @export
// [[Rcpp::export]]
List revert_dist_list_cpp(int n_t, int n_c, std::vector<int> startn,
                                            std::vector<int> endn,
                                            std::vector<double> d) {

  int num_edge = startn.size();
  std::vector<int> start_n_new;
  std::vector<int> end_n_new;
  std::vector<double> d_new;

  for(int i = 0; i < n_c; ++i) {
   std::vector<int> edge_ind;
   std::vector<int> startn_at_edge_ind;
   std::vector<double> d_at_edge_ind;
   for(int j = 0; j < num_edge; ++j) {
     if (endn[j] == i + n_t + 1) {
       edge_ind.push_back(j);
       startn_at_edge_ind.push_back(startn[j] + n_c);
       d_at_edge_ind.push_back(d[j]);
     }
   }

   int num_edge_temp = edge_ind.size();

   int* arr = 0;
   arr = new int[num_edge_temp];

   std::fill_n(arr, num_edge_temp, i + 1);
   start_n_new.insert(start_n_new.end(), arr, arr+num_edge_temp);
   end_n_new.insert(end_n_new.end(), startn_at_edge_ind.begin(), startn_at_edge_ind.end());
   d_new.insert(d_new.end(), d_at_edge_ind.begin(), d_at_edge_ind.end());
  }
  return List::create(start_n_new,end_n_new,d_new);
}




