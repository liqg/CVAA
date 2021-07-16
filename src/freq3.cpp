#include <Rcpp.h>
using namespace Rcpp;


IntegerVector freq3_v (NumericVector x, NumericVector groupEnd, int id, int nr, int nc) {

  int ng = groupEnd.size();
  id = id - 1;

  IntegerVector freq(ng*3*(nc-1), 0);

  int i=0, j=0, k=0, a=0, b=0, c=0, m=0;
  for (i=0; i<nc; i++){
    k=0;
    if(i==id){
      m=1;
      continue;
    }
    for(j=0; j<nr; j++){
      if (x[id*nr+j] < x[i*nr+j] )
        a++;
      else if (x[i*nr+j] == x[id*nr+j])
        b++;
      else
        c++;

      if(j == groupEnd[k]-1){
        freq[(i-m)*ng*3+k*3] = a;
        freq[(i-m)*ng*3+k*3+1] = b;
        freq[(i-m)*ng*3+k*3+2] = c;
        a = b = c = 0;
        k++;
      }
    }
  }
  return (wrap(freq));
}

// [[Rcpp::export]]
IntegerVector freq3 (NumericMatrix x, NumericVector groupEnd, int id) {

  int ng = groupEnd.size();
  int nc = x.ncol();
  int nr = x.nrow();

  id = id - 1; // 0-based

  IntegerVector freq(ng*3*(nc-1), 0);

  int i=0, j=0, k=0, a=0, b=0, c=0, m=0;
  // i: col index, j: row index, k: group index,
  // a, b, c: freq in each group
  // m: after id
  for (i=0; i<nc; i++){
    k=0;
    if(i==id){
      m=1;
      continue;
    }
    for(j=0; j<nr; j++){
      if (x(j, id) < x(j, i))
        a++;
      else if (x(j, i) == x(j, id))
        b++;
      else
        c++;

      if(j == groupEnd[k]-1){ // 0-base
        freq[(i-m)*ng*3+k*3] = a;
        freq[(i-m)*ng*3+k*3+1] = b;
        freq[(i-m)*ng*3+k*3+2] = c;

        a = b = c = 0;
        k++;
      }
    }
  }
  return (wrap(freq));
}

