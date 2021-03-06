#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector c_insert_punctuation_bits(NumericVector punctured_message, NumericVector punctuation_vector, int rows, int cols)
{
	int amount_ones = sum(punctuation_vector);
	int amount_zeros = punctuation_vector.size() - amount_ones;
	int punctured_message_length =  punctured_message.size();

	//Anzahl von 1er bits die beim letzten Matrixdurchlauf verwendet wurden
	int temp = punctured_message_length % amount_ones;

	int i = 0;
	int ones = 0;
	int zeros = 0;
	while(((i % rows) != 0) || ones < temp)
	{
		if(punctuation_vector[i] == 1) {
			ones++;
		} 
		else {
			zeros++;
		}
		i++;
	}

	int additional_bits = (punctured_message_length / amount_ones) * amount_zeros + zeros;

	NumericVector orig_message(punctured_message_length + additional_bits);
	int j = 0;
	for(int i = 0; i < orig_message.size(); i++)
	{
			if(punctuation_vector[i % (rows*cols)] == 1) {
				orig_message[i] = punctured_message[j++];
			} 
			else {
				orig_message[i] = 0;
			}
	}

	return orig_message;
}
