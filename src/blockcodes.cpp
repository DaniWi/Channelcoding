#include <Rcpp.h>
using namespace Rcpp;

#define DEBUG 0


// [[Rcpp::export]]
IntegerVector c_bchEncode
  (	IntegerVector input,
    IntegerVector genPoly,
    int length,
    int k
  ) {


  int    i, j;
  int    feedback;

  int bb[length-k];
  IntegerVector code(length);


  for (i = 0; i < length - k; i++)
    bb[i] = 0;
  for (i = k - 1; i >= 0; i--) {
    feedback = input[i] ^ bb[length - k - 1];
    if (feedback != 0) {
      for (j = length - k - 1; j > 0; j--)
        if (genPoly[j] != 0)
          bb[j] = bb[j - 1] ^ feedback;
        else
          bb[j] = bb[j - 1];
        bb[0] = genPoly[0] && feedback;
    } else {
      for (j = length - k - 1; j > 0; j--)
        bb[j] = bb[j - 1];
      bb[0] = 0;
    }
  }

  for (i = 0; i < length - k; i++)
    code[i] = bb[i];
  for (i = 0; i < k; i++)
    code[i + length - k] = input[i];

  return code;
}

// [[Rcpp::export]]
List c_bchDecode
  (	IntegerVector input,
    int length,
    int m,
    int k,
    int t,
    IntegerVector alpha_to,
    IntegerVector index_of
  ) {

  int i, j, u, q, t2, n, mask, count = 0, syn_error = 0;
  int elp[1026][1024], d[1026], l[1026], u_lu[1026], s[1025];
  int root[200], loc[200], err[1024], reg[201];
  bool failed = false;

  n = 1;
  for (i = 0; i <= m; i++)
    n *= 2;

  n = n / 2 - 1;

  t2 = 2 * t;

  /* first form the syndromes */
  //printf("S(x) = ");

  for (i = 1; i <= t2; i++) {
    s[i] = 0;
    for (j = 0; j < length; j++)
      if (input[j] != 0)
        s[i] ^= alpha_to[(i * j) % n];

      if (s[i] != 0)
        syn_error = 1; /* set error flag if non-zero syndrome */
  /*
      * Note:    If the code is used only for ERROR DETECTION, then
      *          exit program here indicating the presence of errors.
      */
  /* convert syndrome from polynomial form to index form  */
  s[i] = index_of[s[i]];
  //printf("%3d ", s[i]);
  }
  // printf("\n");

  if (syn_error) { /* if there are errors, try to correct them */
  /*
    * Compute the error location polynomial via the Berlekamp
    * iterative algorithm. Following the terminology of Lin and
    * Costello's book :   d[u] is the 'mu'th discrepancy, where
    * u='mu'+1 and 'mu' (the Greek letter!) is the step number
    * ranging from -1 to 2*t (see L&C),  l[u] is the degree of
    * the elp at that step, and u_l[u] is the difference between
    * the step number and the degree of the elp.
    */



  /* initialise table entries */
  d[0] = 0; /* index form */
  d[1] = s[1]; /* index form */
  elp[0][0] = 0; /* index form */
  elp[1][0] = 1; /* polynomial form */
  for (i = 1; i < t2; i++) {
    elp[0][i] = -1; /* index form */
  elp[1][i] = 0; /* polynomial form */
  }
  l[0] = 0;
  l[1] = 0;
  u_lu[0] = -1;
  u_lu[1] = 0;
  u = 0;

  do {
    u++;
    if (d[u] == -1) {
      l[u + 1] = l[u];
      for (i = 0; i <= l[u]; i++) {
        elp[u + 1][i] = elp[u][i];
        elp[u][i] = index_of[elp[u][i]];
      }
    }
    else
      /*
      * search for words with greatest u_lu[q] for
      * which d[q]!=0
      */
    {
      q = u - 1;
      while ((d[q] == -1) && (q > 0))
        q--;
      /* have found first non-zero d[q]  */
      if (q > 0) {
        j = q;
        do {
          j--;
          if ((d[j] != -1) && (u_lu[q] < u_lu[j]))
            q = j;
        } while (j > 0);
      }

      /*
      * have now found q such that d[u]!=0 and
      * u_lu[q] is maximum
      */
      /* store degree of new elp polynomial */
      if (l[u] > l[q] + u - q)
        l[u + 1] = l[u];
      else
        l[u + 1] = l[q] + u - q;

      /* form new elp(x) */
      for (i = 0; i < t2; i++)
        elp[u + 1][i] = 0;
      for (i = 0; i <= l[q]; i++)
        if (elp[q][i] != -1)
          elp[u + 1][i + u - q] = alpha_to[(d[u] + n - d[q] + elp[q][i]) % n];
        for (i = 0; i <= l[u]; i++) {
          elp[u + 1][i] ^= elp[u][i];
          elp[u][i] = index_of[elp[u][i]];
        }
    }
    u_lu[u + 1] = u - l[u + 1];

    /* form (u+1)th discrepancy */
    if (u < t2) {
      /* no discrepancy computed on last iteration */
      if (s[u + 1] != -1)
        d[u + 1] = alpha_to[s[u + 1]];
      else
        d[u + 1] = 0;
      for (i = 1; i <= l[u + 1]; i++)
        if ((s[u + 1 - i] != -1) && (elp[u + 1][i] != 0))
          d[u + 1] ^= alpha_to[(s[u + 1 - i]
                                  + index_of[elp[u + 1][i]])
          % n];
          /* put d[u+1] into index form */
          d[u + 1] = index_of[d[u + 1]];
    }
  } while ((u < t2) && (l[u + 1] <= t));

  u++;
  if (l[u] <= t) { /* Can correct errors */
          /* put elp into index form */
          for (i = 0; i <= l[u]; i++)
            elp[u][i] = index_of[elp[u][i]];

   /* printf("sigma(x) = ");
    for (i = 0; i <= l[u]; i++)
      printf("%3d ", elp[u][i]);
    printf("\n");
    printf("Roots: ");
*/
    /* Chien search: find roots of the error location polynomial */
    for (i = 1; i <= l[u]; i++)
      reg[i] = elp[u][i];
    count = 0;
    for (i = 1; i <= n; i++) {
      q = 1;
      for (j = 1; j <= l[u]; j++)
        if (reg[j] != -1) {
          reg[j] = (reg[j] + j) % n;
          q ^= alpha_to[reg[j]];
        }
        if (!q) { /* store root and error
          * location number indices */
    root[count] = i;
          loc[count] = n - i;
          count++;
          //printf("%3d ", n - i);
        }
    }
    //printf("\n");
    if (count == l[u])
      /* no. roots = degree of elp hence <= t errors */
      for (i = 0; i < l[u]; i++) {
        input[loc[i]] ^= 1;
        loc[i] = 0;
      }
      else /* elp has degree >t hence cannot solve */
        failed = true;
      //printf("Incomplete decoding: errors detected\n");
  }
  }

  IntegerVector decoded(k);
  for(i=0;i<k;i++)
    decoded[i] = input[length-k+i];

  List ret = List::create(Rcpp::Named("decoded") = decoded,
                          Rcpp::Named("input.corrected") = input,
                          Rcpp::Named("failed") = failed);
  return ret;



}

// [[Rcpp::export]]
List c_getGeneratorPoly
  (	int length,
    int m,
    int t
  ) {

  int n, k, d;
  int p[m];
  int ii, jj, ll, kaux, i, mask, ninf;
  int test, aux, nocycles, root, noterms, rdncy;
  int cycle[1024][m], size[1024], min[1024], zeros[1024];

  IntegerVector alpha_to(pow(2,m));
  IntegerVector index_of(pow(2,m));

  IntegerVector g(length);

  for (i = 1; i < m; i++)
    p[i] = 0;
  p[0] = p[m] = 1;
  if (m == 2)
    p[1] = 1;
  else if (m == 3)
    p[1] = 1;
  else if (m == 4)
    p[1] = 1;
  else if (m == 5)
    p[2] = 1;
  else if (m == 6)
    p[1] = 1;
  else if (m == 7)
    p[1] = 1;
  else if (m == 8)
    p[4] = p[5] = p[6] = 1;
  else if (m == 9)
    p[4] = 1;
  else if (m == 10)
    p[3] = 1;
  else if (m == 11)
    p[2] = 1;
  else if (m == 12)
    p[3] = p[4] = p[7] = 1;
  else if (m == 13)
    p[1] = p[3] = p[4] = 1;
  else if (m == 14)
    p[1] = p[11] = p[12] = 1;
  else if (m == 15)
    p[1] = 1;
  else if (m == 16)
    p[2] = p[3] = p[5] = 1;
  else if (m == 17)
    p[3] = 1;
  else if (m == 18)
    p[7] = 1;
  else if (m == 19)
    p[1] = p[5] = p[6] = 1;
  else if (m == 20)
    p[3] = 1;

  n = 1;
  for (i = 0; i <= m; i++)
    n *= 2;

  n = n / 2 - 1;
  ninf = (n + 1) / 2 - 1;

  //Construct GF(2^m)
  mask = 1;
  alpha_to[m] = 0;
  for (i = 0; i < m; i++) {
    alpha_to[i] = mask;
    index_of[alpha_to[i]] = i;
    if (p[i] != 0)
      alpha_to[m] ^= mask;
    mask <<= 1;
  }
  index_of[alpha_to[m]] = m;
  mask >>= 1;
  for (i = m + 1; i < n; i++) {
    if (alpha_to[i - 1] >= mask)
      alpha_to[i] = alpha_to[m] ^ ((alpha_to[i - 1] ^ mask) << 1);
    else
      alpha_to[i] = alpha_to[i - 1] << 1;
    index_of[alpha_to[i]] = i;
  }
  index_of[0] = -1;

  /* Generate cycle sets modulo n, n = 2**m - 1 */
  cycle[0][0] = 0;
  size[0] = 1;
  cycle[1][0] = 1;
  size[1] = 1;
  jj = 1; /* cycle set index */

  do {
    /* Generate the jj-th cycle set */
    ii = 0;
    do {
      ii++;
      cycle[jj][ii] = (cycle[jj][ii - 1] * 2) % n;
      size[jj]++;
      aux = (cycle[jj][ii] * 2) % n;
    } while (aux != cycle[jj][0]);
    /* Next cycle set representative */
    ll = 0;
    do {
      ll++;
      test = 0;
      for (ii = 1; ((ii <= jj) && (!test)); ii++)
        /* Examine previous cycle sets */
        for (kaux = 0; ((kaux < size[ii]) && (!test)); kaux++)
          if (ll == cycle[ii][kaux])
            test = 1;
    } while ((test) && (ll < (n - 1)));
    if (!(test)) {
      jj++; /* next cycle set index */
        cycle[jj][0] = ll;
        size[jj] = 1;
    }
  } while (ll < (n - 1));
  nocycles = jj; /* number of cycle sets modulo n */

        d = 2 * t + 1;

        /* Search for roots 1, 2, ..., d-1 in cycle sets */
        kaux = 0;
        rdncy = 0;
        for (ii = 1; ii <= nocycles; ii++) {
          min[kaux] = 0;
          test = 0;
          for (jj = 0; ((jj < size[ii]) && (!test)); jj++)
            for (root = 1; ((root < d) && (!test)); root++)
              if (root == cycle[ii][jj]) {
                test = 1;
                min[kaux] = ii;
              }
              if (min[kaux]) {
                rdncy += size[min[kaux]];
                kaux++;
              }
        }
        noterms = kaux;
        kaux = 1;
        for (ii = 0; ii < noterms; ii++)
          for (jj = 0; jj < size[min[ii]]; jj++) {
            zeros[kaux] = cycle[min[ii]][jj];
            kaux++;
          }

          k = length - rdncy;

        if (k < 0) {
          Rprintf("Parameters invalid!\n");
        }

        //  printf("This is a (%d, %d, %d) binary BCH code\n", length, k, d);


        /* Compute the generator polynomial */
        g[0] = alpha_to[zeros[1]];
        g[1] = 1; /* g(x) = (X + zeros[1]) initially */
        for (ii = 2; ii <= rdncy; ii++) {
          g[ii] = 1;
          for (jj = ii - 1; jj > 0; jj--)
            if (g[jj] != 0)
              g[jj] = g[jj - 1] ^ alpha_to[(index_of[g[jj]] + zeros[ii]) % n];
            else
              g[jj] = g[jj - 1];
            g[0] = alpha_to[(index_of[g[0]] + zeros[ii]) % n];
        }

        IntegerVector g_cleaned(length-k+1);
        for(i=0;i<=length-k;i++)
          g_cleaned[i] = g[i];

        List ret = List::create(Rcpp::Named("gen.poly") = g_cleaned,
                                Rcpp::Named("data.length") = k,
                                Rcpp::Named("alpha_to") = alpha_to,
                                Rcpp::Named("index_of") = index_of);
        return ret;
}

