#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rcpp.h>
using namespace Rcpp;

#define DEBUG 0
//noise standard deviation
double sigma = 1.0;
double Lc = 2/(sigma * sigma);

NumericVector c_sova
(
	NumericVector x_d,
	NumericVector p_d,
	NumericVector La,
	int	is_terminated,
	int N,
	int M,
	IntegerMatrix previous_state,
	IntegerMatrix output,
	int output_index
)
{
	const int msg_len = x_d.size();		// includes termination bits (M termination bits)
	const int NUM_STATES = pow(2,M);

	double deltamin;
	double delta[msg_len+1][NUM_STATES];
	double metric[msg_len+1][NUM_STATES];
	int survivor_bit[msg_len+1][NUM_STATES];
	int survivor_states[msg_len+1];
	int previous_matrix_column[msg_len+1][NUM_STATES];
	// NOTE: arrays have dimension msg_len+M+1 because i.e. delta[1][state] is delta from symbol at t=0 to t=1
	//		 so delta from symbol at t=msg_len+M-1 to t=msg_len+M (+M because of termination) is stored in
	//		 delta[msg_len+M][state] -> therefore delta[msg_len+M+1][NUM_STATES]
	//		 index 0 is not used!

	// initialization
	metric[0][0] = 0;
	for (int i = 1; i < NUM_STATES; i++)
	{
		metric[0][i] = -1000;
	}
	for (int i = 0; i < msg_len+1; i++)
	{
		for (int j = 0; j < NUM_STATES; j++)
		{
			previous_matrix_column[i][j] = -1;
		}
	}

	// loop: time t
	for (int t = 1; t < msg_len+1; t++)
	{
		// loop: state s
		for (int s = 0; s < NUM_STATES; s++)
		{
			double Max[2];

			if (previous_state(s, 0) == -1 || previous_state(s, 1) == -1)
			{
				// Max[0] contains the metrik for previous State from column 0 and 1
				// Max[1] contains the metrik for previous State from column 2
				if (previous_state(s, 0) == -1)
				{
					// no previous for input bit 0
					survivor_bit[t][s] = 1;
				}
				else if (previous_state(s, 1) == -1)
				{
					// no previous for input bit 1
					survivor_bit[t][s] = 0;
				}
				Max[0] = metric[t-1][previous_state(s,survivor_bit[t][s])];
				Max[1] = metric[t-1][previous_state(s,2)];

				//calculate the first possibility
				Max[0] += (1 - 2 * ((output(previous_state(s,survivor_bit[t][s]),survivor_bit[t][s]) >> (N-1)) & 0x01)) * (Lc * x_d[t-1] + La[t-1]);
				Max[0] += (1 - 2 * ((output(previous_state(s,survivor_bit[t][s]),survivor_bit[t][s]) >> (N-output_index)) & 0x01)) * Lc * p_d[t-1];

				//calculate the second possibility
				Max[1] += (1 - 2 * ((output(previous_state(s,2),survivor_bit[t][s]) >> (N-1)) & 0x01)) * (Lc * x_d[t-1] + La[t-1]);
				Max[1] += (1 - 2 * ((output(previous_state(s,2),survivor_bit[t][s]) >> (N-output_index)) & 0x01)) * Lc * p_d[t-1];

				previous_matrix_column[t][s] = (Max[0] > Max[1]) ? survivor_bit[t][s] : 2;
			}
			else
			{
				// there are previous states for input bits 0 and 1
				// Max[0] contains the metrik for previous State through previous bit 0
				// Max[1] contains the metrik for previous State through previous bit 1
				for (int i = 0; i < 2; i++)
				{
					Max[i] = metric[t-1][previous_state(s,i)];

					Max[i] += (1 - 2 * ((output(previous_state(s,i),i) >> (N-1)) & 0x01)) * (Lc * x_d[t-1] + La[t-1]);
					Max[i] += (1 - 2 * ((output(previous_state(s,i),i) >> (N-output_index)) & 0x01)) * Lc * p_d[t-1];
				}
				survivor_bit[t][s] = (Max[0] > Max[1]) ? 0 : 1;
			}

			if (Max[0] > Max[1])
			{
				metric[t][s] = Max[0];
				delta[t][s] = (Max[0] - Max[1]) / 2;
			}
			else
			{
				metric[t][s] = Max[1];
				delta[t][s] = (Max[1] - Max[0]) / 2;
			}

			#if DEBUG > 1
			Rprintf("\n");
			Rprintf("Metric[%d][%d] = % 3.2f\t delta[%d][%d] = % 3.2f\t M0 = % 3.2f\t M1 = % 3.2f\t survivor_bit = %d", t, s, metric[t][s], t, s, delta[t][s], Max[0], Max[1], survivor_bit[t][s]);
			#endif
		}
	}

	if (is_terminated > 0)
	{
		survivor_states[msg_len] = 0;
	}
	else
	{
		double max_metric = metric[msg_len][0];
		survivor_states[msg_len] = 0;
		for (int s = 1; s < NUM_STATES; s++)
		{
			if (max_metric < metric[msg_len][s])
			{
				max_metric = metric[msg_len][s];
				survivor_states[msg_len] = s;
			}
		}
	}

	// reconstruction of survivor states, from right to left
	// survivor state at index msg_len is known, start loop at 2nd but last index (= msg_len-1)
	for (int t = msg_len-1; t >= 0; t--)
	{
		if (previous_matrix_column[t+1][survivor_states[t+1]] == -1)
		{
			survivor_states[t] = previous_state(survivor_states[t+1], survivor_bit[t+1][survivor_states[t+1]]);
		}
		else
		{
			survivor_states[t] = previous_state(survivor_states[t+1], previous_matrix_column[t+1][survivor_states[t+1]]);
		}

		#if DEBUG > 1
		Rprintf("\n");
		Rprintf("survivor_states[%d] = %d", t, survivor_states[t]);
		#endif
	}



	int s;
	for (int t = M+1; t < msg_len+1; t++)
	{
		deltamin = delta[t][survivor_states[t]];
		// s ... state which is the second possibility to come to the survivor state

		if (previous_matrix_column[t][survivor_states[t]] == -1)
		{
			// if no decision was made, bit has to flip
			s = previous_state(survivor_states[t], survivor_bit[t][survivor_states[t]] * (-1) + 1);
		}
		else
		{
			// if decision was made, change column
			// 2 cases:	a) 0 <--> 2
			//			b) 1 <--> 2
			// 0 or 1 comes from the survivor bit
			int column;
			if (previous_matrix_column[t][survivor_states[t]] == 2)
			{
				column = survivor_bit[t][survivor_states[t]];
			}
			else
			{
				column = 2;
			}
			s = previous_state(survivor_states[t], column);
		}

		for (int i = t-1; i > 0; i--)
		{
			//searches the smallest delta from survivor state
			if (delta[i][survivor_states[i]] < deltamin)
			{
				deltamin = delta[i][survivor_states[i]];
			}
			if (survivor_bit[i][survivor_states[i]] != survivor_bit[i][s])
			{
				delta[i][survivor_states[i]] = deltamin;
				#if DEBUG == 3
				Rprintf("\n");
				Rprintf("update\t t  = %3d\t i = %3d\t delta = %3.2f", t, i, delta[i][survivor_states[i]]);
				#endif
			}

			if (previous_matrix_column[i][s] == -1)
			{
				//if no decision was made, bit has to flip
				s = previous_state(s, survivor_bit[i][s]);
			}
			else
			{
				//if decision was made, switch column
				s = previous_state(s, previous_matrix_column[i][s]);
			}
		}

	}

	NumericVector soft_output(msg_len);

	for(int t = 1; t < msg_len+1; t++)
	{
		soft_output[t-1] = delta[t][survivor_states[t]] * (1 - 2 * survivor_bit[t][survivor_states[t]]) - La[t-1] - Lc * x_d[t-1];
	}

	return soft_output;
}

// [[Rcpp::export]]
List c_turbo_decode
(
	NumericVector x_noisy,
	NumericVector parity_noisy1,
	NumericVector parity_noisy2,
	IntegerVector permutation,
	int N_ITERATION,
	int N,
	int M,
	IntegerMatrix previous_state,
	IntegerMatrix output,
	int output_index
)
{
	const int msg_len = x_noisy.size();

	//lists to save the temporary results during soft decoding
	List decode1(N_ITERATION);
	List decode1I(N_ITERATION);
	List decode2Back(N_ITERATION);
	List decode2IBack(N_ITERATION);
	List tempResultSoft(N_ITERATION);
	List tempResultHard(N_ITERATION);

	NumericVector x_d_p(msg_len);  	//noisy data permutated
	NumericVector Le1;    			//decoder #1 extrinsic likelihood
	NumericVector Le1_p(msg_len);  	//decoder #1 extrinsic likelihood permuted
	NumericVector Le2;    			//decoder #2 extrinsic likelihood
	NumericVector Le2_ip(msg_len); 	//decoder #2 extrinsic likelihood inverse permuted


    for(int k = 0; k < msg_len; k++)
	{
		Le2_ip[k] = 0;
	}

    for(int i = 0; i < N_ITERATION; i++)
    {
		Le1 = c_sova(x_noisy, parity_noisy1, Le2_ip, 1, N, M, previous_state, output, output_index);

       //permute decoder#1 likelihoods to match decoder#2 order
    	for(int k = 0; k < msg_len; k++)
    	{
			Le1_p[k] = Le1[permutation[k]];
			x_d_p[k] = x_noisy[permutation[k]];
    	}

		Le2 = c_sova(x_d_p, parity_noisy2, Le1_p, 0, N, M, previous_state, output, output_index);

        //inverse permute decoder#2 likelihoods to match decoder#1 order
    	for(int k = 0; k < msg_len; k++)
		{
    		Le2_ip[permutation[k]] = Le2[k];
		}

    	#if DEBUG > 0
    	Rprintf("\nIteration %d:\n",i+1);
		for(int k = 0; k < msg_len; k++)
		{
 			Rprintf("Le1[%i] = % 3.2f\t", k, Le1[k]);
 			Rprintf("Le1_p[%i] = % 3.2f\t",k, Le1_p[k]);
 			Rprintf("Le2[%i] = % 3.2f\t",k, Le2[k]);
 			Rprintf("Le2_ip[%i] = % 3.2f\n", k, Le2_ip[k]);
		}
		Rprintf("\n");
		#endif

		//vectors to save in the list
		NumericVector Le1_o(msg_len);
		NumericVector Le1_p_o(msg_len);
		NumericVector Le2_o(msg_len);
		NumericVector Le2_ip_o(msg_len);
		NumericVector soft(msg_len);
		NumericVector hard(msg_len);

		//copy the vectors
		for(int k = 0; k < msg_len; k++)
		{
 			Le1_o[k] = Le1[k];
			Le1_p_o[k] = Le1_p[k];
			Le2_o[k] = Le2[k];
			Le2_ip_o[k] = Le2_ip[k];
			soft[k] = Lc * x_noisy[k] + Le1[k] + Le2_ip[k];
			hard[k] = (soft[k] >= 0.0) ? 0 : 1;
		}

		//save the vectors into the list
		decode1[i] = Le1_o;
		decode1I[i] = Le1_p_o;
		decode2Back[i] = Le2_o;
		decode2IBack[i] = Le2_ip_o;
		tempResultSoft[i] = soft;
		tempResultHard[i] = hard;

	}

	NumericVector soft_output(msg_len);
	IntegerVector hard_output(msg_len);

    //calculate overall likelihoods and then calculate hard decisions
    for(int k = 0; k < msg_len; k++)
    {
    	soft_output[k] = Lc * x_noisy[k] + Le1[k] + Le2_ip[k]; 		//soft decision
    	hard_output[k] = (soft_output[k] >= 0.0) ? 0 : 1;         //hard decision
	}

	#if DEBUG > 0
    //print soft decisions
    Rprintf("\nSoftDecision\n");
	for(int k = 0; k < msg_len ; k++)
		Rprintf("L_h[%i] = % .2f\n", k, soft_output[k]);
	Rprintf("\n");

    //print hard decisions
	Rprintf("output_hard = ");
    for(int k = 0; k < msg_len; k++)
    	Rprintf("%i", hard_output[k]);
    Rprintf("\n");
    #endif

    List infoDisp = List::create(Rcpp::Named("origI") = x_d_p ,
								 Rcpp::Named("decode1") = decode1,
								 Rcpp::Named("decode1I") = decode1I,
								 Rcpp::Named("decode2Back") = decode2Back,
								 Rcpp::Named("decode2IBack") = decode2IBack,
								 Rcpp::Named("tempResultSoft") = tempResultSoft,
								 Rcpp::Named("tempResultHard") = tempResultHard);


    List result = List::create(Rcpp::Named("soft.output") = soft_output,
							   Rcpp::Named("hard.output") = hard_output,
							   Rcpp::Named("disp.info") = infoDisp);

	return result;
}
