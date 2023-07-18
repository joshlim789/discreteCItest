#include "pDCorTest.h"

// [[Rcpp::export]]
StringVector matrix_to_string (StringMatrix sep_vectors)
{
  int row_n = sep_vectors.nrow();
  StringVector con_string (row_n);
  for (int i = 0; i < row_n; i++)
  {
    //Concatenates each row as a singular element of a vector
    con_string[i] = collapse(sep_vectors(i,_));
  }
  return con_string;
}

// [[Rcpp::export]]
double get_G2_one(StringVector A, StringVector B, int tot_Au_size, int tot_Bu_size)
{
  //Contingency Table
  StringVector A_uniq = unique(A);
  StringVector B_uniq = unique(B);
  
  int Au_size = A_uniq.size();
  int Bu_size = B_uniq.size();
  
  int A_size = A.size();
  
  //Counting each combination
  arma::mat O (tot_Au_size, tot_Bu_size);
  
  for (int m = 0; m < Au_size; m ++)
  {
    for (int n = 0; n < Bu_size; n ++)
    {
      for (int i = 0; i < A_size; i ++)
      {
        if((A[i] == A_uniq[m]) & (B[i] == B_uniq[n]))
        {
          O (m, n) ++;
        }
      }
    }
  }

  
  //Expected counts
  arma::rowvec rowz = arma::sum(O, 0);
  arma::colvec colz= arma::sum(O, 1);
  arma::mat E =  colz * rowz / arma::accu(O);
  
  
  //G2 Adding, accounting for 0 observed counts
  double G2 = 0;
  
  for(int m = 0; m < tot_Au_size; m ++)
  {
    for(int n = 0; n < tot_Bu_size; n ++)
    {
      if(O(m, n) != 0)
      {
        G2 = G2 + O(m, n)*(log(O(m, n)/E(m, n))); 
      }
    }
  }
  
  return 2*G2;
}

// [[Rcpp::export]]
double get_G2_all(StringVector A, StringVector B, StringVector S)
{
  
  StringVector S_uniq = unique(S);
  int Su_size = S_uniq.size();
  arma::vec G_squares (Su_size);
  
  int tot_Au_size = unique(A).size();
  int tot_Bu_size = unique(B).size();
  
  //Calculate for each level of S
  for (int i = 0; i < Su_size; i++)
  {
    
    StringVector A_sub;
    StringVector B_sub;
    
    for (int j = 0; j < S.size(); j++)
    {
      if (S[j] == S_uniq[i])
      {
        A_sub.push_back(A[j]);
        B_sub.push_back(B[j]);
      }
    }
    G_squares[i] = get_G2_one(A_sub, B_sub, tot_Au_size, tot_Bu_size);
  }
  
  double G_stat = arma::sum(G_squares);
  return G_stat;
}

// [[Rcpp::export]]
List condInttestdis(StringMatrix df, const size_t &i,const size_t &j,
                     const arma::uvec &k, const double &signif_level)
{
  size_t k_size = k.size();
  
  //Setting up vectors and conditioning set
  StringVector A = df(_,i);
  StringVector B = df(_,j);
  StringMatrix S_m(df.nrow(), k_size);
  
  int S_df = 1;
  
  //Looping through to get the elements and df of conditioning set
  for (size_t j = 0; j < k_size; j++)
  {
    StringVector con_col = df(_ , k(j));
    S_m(_, j) = con_col;
    S_df = S_df * unique(con_col).size();
  }
  
  //Convert to one vector to examine the different vector levels
  StringVector S = matrix_to_string(S_m);
  
  //Calculating the relevant info
  double statistic = get_G2_all(A, B, S);
  int dof = (unique(A).size() - 1) * (unique(B).size() - 1) * S_df;
  double cutoff = R::qchisq(1-signif_level,dof,true,false);
  bool accept_H0 = std::abs(statistic) <= cutoff;
  double pval = 1 - R::pchisq(statistic, dof, true,false);
  return List::create(
    _["result"]=accept_H0,
    _["statistic"]=statistic,
    _["pval"]=pval
  );
}
