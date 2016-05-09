// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// c_bchEncode
IntegerVector c_bchEncode(IntegerVector input, IntegerVector genPoly, int length, int k);
RcppExport SEXP channelcoding_c_bchEncode(SEXP inputSEXP, SEXP genPolySEXP, SEXP lengthSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type input(inputSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type genPoly(genPolySEXP);
    Rcpp::traits::input_parameter< int >::type length(lengthSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    __result = Rcpp::wrap(c_bchEncode(input, genPoly, length, k));
    return __result;
END_RCPP
}
// c_bchDecode
IntegerVector c_bchDecode(IntegerVector input, IntegerVector alpha_to, IntegerVector index_of, int length, int m, int t, int k);
RcppExport SEXP channelcoding_c_bchDecode(SEXP inputSEXP, SEXP alpha_toSEXP, SEXP index_ofSEXP, SEXP lengthSEXP, SEXP mSEXP, SEXP tSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type input(inputSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type alpha_to(alpha_toSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type index_of(index_ofSEXP);
    Rcpp::traits::input_parameter< int >::type length(lengthSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    __result = Rcpp::wrap(c_bchDecode(input, alpha_to, index_of, length, m, t, k));
    return __result;
END_RCPP
}
// c_getGeneratorPoly
List c_getGeneratorPoly(int length, int m, int t);
RcppExport SEXP channelcoding_c_getGeneratorPoly(SEXP lengthSEXP, SEXP mSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type length(lengthSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    __result = Rcpp::wrap(c_getGeneratorPoly(length, m, t));
    return __result;
END_RCPP
}
// c_generateMatrices
List c_generateMatrices(int N, int M, IntegerVector generator);
RcppExport SEXP channelcoding_c_generateMatrices(SEXP NSEXP, SEXP MSEXP, SEXP generatorSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type generator(generatorSEXP);
    __result = Rcpp::wrap(c_generateMatrices(N, M, generator));
    return __result;
END_RCPP
}
// c_generateMatrices_rsc
List c_generateMatrices_rsc(int N, int M, IntegerVector generator);
RcppExport SEXP channelcoding_c_generateMatrices_rsc(SEXP NSEXP, SEXP MSEXP, SEXP generatorSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type generator(generatorSEXP);
    __result = Rcpp::wrap(c_generateMatrices_rsc(N, M, generator));
    return __result;
END_RCPP
}
// c_convolutionEncode
IntegerVector c_convolutionEncode(IntegerVector input, int N, int M, IntegerMatrix nextState, IntegerMatrix output, int rsc, IntegerVector termination, int terminate);
RcppExport SEXP channelcoding_c_convolutionEncode(SEXP inputSEXP, SEXP NSEXP, SEXP MSEXP, SEXP nextStateSEXP, SEXP outputSEXP, SEXP rscSEXP, SEXP terminationSEXP, SEXP terminateSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type input(inputSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type nextState(nextStateSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type output(outputSEXP);
    Rcpp::traits::input_parameter< int >::type rsc(rscSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type termination(terminationSEXP);
    Rcpp::traits::input_parameter< int >::type terminate(terminateSEXP);
    __result = Rcpp::wrap(c_convolutionEncode(input, N, M, nextState, output, rsc, termination, terminate));
    return __result;
END_RCPP
}
// c_convolutionDecode
List c_convolutionDecode(NumericVector code, int N, int M, IntegerMatrix previousState, IntegerMatrix output, int IsTerminated);
RcppExport SEXP channelcoding_c_convolutionDecode(SEXP codeSEXP, SEXP NSEXP, SEXP MSEXP, SEXP previousStateSEXP, SEXP outputSEXP, SEXP IsTerminatedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type code(codeSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type previousState(previousStateSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type output(outputSEXP);
    Rcpp::traits::input_parameter< int >::type IsTerminated(IsTerminatedSEXP);
    __result = Rcpp::wrap(c_convolutionDecode(code, N, M, previousState, output, IsTerminated));
    return __result;
END_RCPP
}
// c_convolutionDecode_hard
List c_convolutionDecode_hard(IntegerVector code, int N, int M, IntegerMatrix previousState, IntegerMatrix output, int IsTerminated);
RcppExport SEXP channelcoding_c_convolutionDecode_hard(SEXP codeSEXP, SEXP NSEXP, SEXP MSEXP, SEXP previousStateSEXP, SEXP outputSEXP, SEXP IsTerminatedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type code(codeSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type previousState(previousStateSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type output(outputSEXP);
    Rcpp::traits::input_parameter< int >::type IsTerminated(IsTerminatedSEXP);
    __result = Rcpp::wrap(c_convolutionDecode_hard(code, N, M, previousState, output, IsTerminated));
    return __result;
END_RCPP
}
// gcd_polynomial
int gcd_polynomial(IntegerVector x);
RcppExport SEXP channelcoding_gcd_polynomial(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    __result = Rcpp::wrap(gcd_polynomial(x));
    return __result;
END_RCPP
}
// c_insert_punctuation_bits
NumericVector c_insert_punctuation_bits(NumericVector punctured_message, NumericVector punctuation_vector, int rows, int cols);
RcppExport SEXP channelcoding_c_insert_punctuation_bits(SEXP punctured_messageSEXP, SEXP punctuation_vectorSEXP, SEXP rowsSEXP, SEXP colsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type punctured_message(punctured_messageSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type punctuation_vector(punctuation_vectorSEXP);
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    __result = Rcpp::wrap(c_insert_punctuation_bits(punctured_message, punctuation_vector, rows, cols));
    return __result;
END_RCPP
}
// c_turbo_decode
List c_turbo_decode(NumericVector x_noisy, NumericVector parity_noisy1, NumericVector parity_noisy2, IntegerVector permutation, int N_ITERATION, int N, int M, IntegerMatrix previous_state, IntegerMatrix output, int output_index);
RcppExport SEXP channelcoding_c_turbo_decode(SEXP x_noisySEXP, SEXP parity_noisy1SEXP, SEXP parity_noisy2SEXP, SEXP permutationSEXP, SEXP N_ITERATIONSEXP, SEXP NSEXP, SEXP MSEXP, SEXP previous_stateSEXP, SEXP outputSEXP, SEXP output_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x_noisy(x_noisySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type parity_noisy1(parity_noisy1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type parity_noisy2(parity_noisy2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type permutation(permutationSEXP);
    Rcpp::traits::input_parameter< int >::type N_ITERATION(N_ITERATIONSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type previous_state(previous_stateSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type output(outputSEXP);
    Rcpp::traits::input_parameter< int >::type output_index(output_indexSEXP);
    __result = Rcpp::wrap(c_turbo_decode(x_noisy, parity_noisy1, parity_noisy2, permutation, N_ITERATION, N, M, previous_state, output, output_index));
    return __result;
END_RCPP
}
