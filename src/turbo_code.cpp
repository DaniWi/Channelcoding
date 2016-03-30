#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rcpp.h>
using namespace Rcpp;

#define MAP_SOVA 1		//set which algorithm to use [0=>MAP,1=>SOVA]
#define DEBUG 1
//noise standard deviation
double sigma = 1.0;	//Lc=2/(sigma * sigma)=2
double Lc = 2/(sigma * sigma);



/* c_convolutionDecode
 *
 * decodes a code and returns the decoded message (List of soft and hard values)
 * soft decision decoding
 * metric: scalar product (soft value)
 *
 * params
 * code: the code to be decoded (soft values)
 * N: number of output symbols per input symbol
 * M: constraint length, number of memory elements
 * previousState: the previousState-matrix of the convolutional encoder
 * output: the output-matrix of the convolutional encoder
 */

NumericVector c_sova
(
	NumericVector x_d,
	NumericVector p_d,
	NumericVector La,
	int	IsTerminated,
	int N,
	int M,
	IntegerMatrix previousState,
	IntegerMatrix output,
	int output_index
)
{

	const int msgLen = x_d.size();		// includes termination bits (M termination bits)
	const int NUM_STATES = pow(2,M);

	double deltamin;
	double delta[msgLen+1][NUM_STATES];
	double metric[msgLen+1][NUM_STATES];
	int survivorBit[msgLen+1][NUM_STATES];
	int survivorStates[msgLen+1];
	int previousMatrixColumn[msgLen+1][NUM_STATES];
	// NOTE: arrays have dimension msgLen+M+1 because i.e. delta[1][state] is delta from symbol at t=0 to t=1
	//		 so delta from symbol at t=msgLen+M-1 to t=msgLen+M (+M because of termination) is stored in
	//		 delta[msgLen+M][state] -> therefore delta[msgLen+M+1][NUM_STATES]
	//		 index 0 is not used!

	// initialization
	metric[0][0] = 0;
	for (int i = 1; i < NUM_STATES; i++)
	{
		metric[0][i] = -1000;
	}
	for (int i = 0; i < msgLen+1; i++)
	{
		for (int j = 0; j < NUM_STATES; j++)
		{
			previousMatrixColumn[i][j] = -1;
		}
	}

	int Lc = 1;	// ACHTUNG!!!!!!
	int index = 0;	// index to select correct code bit

	// loop: time t
	for (int t = 1; t < msgLen+1; t++) {
		// loop: state s
		for (int s = 0; s < NUM_STATES; s++)
		{
			double Max[2];

			if (previousState(s,0) == -1 || previousState(s,1) == -1)
			{
				// Max[0] erh?lt die Metrik f?r previous State aus Spalte 0 bzw 1
				// Max[1] erh?lt die Metrik f?r previous State aus Spalte 2
				if (previousState(s,0) == -1)
				{
					// no previous for input bit 0
					survivorBit[t][s] = 1;
				}
				else if (previousState(s,1) == -1)
				{
					// no previous for input bit 1
					survivorBit[t][s] = 0;
				}
				Max[0] = metric[t-1][previousState(s,survivorBit[t][s])];
				Max[1] = metric[t-1][previousState(s,2)];

				//Berechnung von erster M?glichkeit
				Max[0] += (((output(previousState(s,survivorBit[t][s]),survivorBit[t][s]) >> (N-1)) & 0x01) * 2 - 1) * (Lc * x_d[t-1] + La[t-1]);
				Max[0] += (((output(previousState(s,survivorBit[t][s]),survivorBit[t][s]) >> (N-output_index)) & 0x01) * 2 - 1) * Lc * p_d[t-1];

				//Berechnung von zweiter M?glichkeit
				Max[1] += (((output(previousState(s,2),survivorBit[t][s]) >> (N-1)) & 0x01) * 2 - 1) * (Lc * x_d[t-1] + La[t-1]);
				Max[1] += (((output(previousState(s,2),survivorBit[t][s]) >> (N-output_index)) & 0x01) * 2 - 1) * Lc * p_d[t-1];

				previousMatrixColumn[t][s] = (Max[0] > Max[1]) ? survivorBit[t][s] : 2;
			}
			else
			{
				// there are previous states for input bits 0 and 1
				// Max[0] erh?lt die Metrik f?r previous State mittels previous bit 0
				// Max[1] erh?lt die Metrik f?r previous State mittels previous bit 1
				for (int i = 0; i < 2; i++)
				{
					Max[i] = metric[t-1][previousState(s,i)];

					Max[i] += (((output(previousState(s,i),i) >> (N-1)) & 0x01) * 2 - 1) * (Lc * x_d[t-1] + La[t-1]);
					Max[i] += (((output(previousState(s,i),i) >> (N-output_index)) & 0x01) * 2 - 1) * Lc * p_d[t-1];
				}

				survivorBit[t][s] = (Max[0] > Max[1]) ? 0 : 1;
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

			#if DEBUG == 1
			printf("\ndelta[%d][%d]=%f\t M0=%.2f \tM1=%.2f \tmax=%f",t,s,delta[t][s], Max[0],Max[1], metric[t][s]);
			#endif
		}

		index += N;
	}

	if (IsTerminated > 0)
	{
		survivorStates[msgLen] = 0;
	}
	else
	{
		double max_metric = metric[msgLen][0];
		survivorStates[msgLen] = 0;
		for (int s = 1; s < NUM_STATES; s++)
		{
			if (max_metric < metric[msgLen][s])
			{
				max_metric = metric[msgLen][s];
				survivorStates[msgLen] = s;
			}
		}
	}

	// reconstruction of survivor states, from right to left
	// survivor state at index msgLen is known, start loop at 2nd but last index (= msgLen-1)
	for (int t = msgLen-1; t >= 0; t--)
	{
		if (previousMatrixColumn[t+1][survivorStates[t+1]] == -1)
		{
			survivorStates[t] = previousState(survivorStates[t+1], survivorBit[t+1][survivorStates[t+1]]);
		}
		else
		{
			//survivorStates[t] = previousState[previousMatrixColumn[t+1][survivorStates[t+1]]][survivorStates[t+1]][survivorBit[t+1][survivorStates[t+1]]];
			survivorStates[t] = previousState(survivorStates[t+1], previousMatrixColumn[t+1][survivorStates[t+1]]);
		}
	}

	int s;
	for (int t = M+1; t < msgLen+1; t++)
	{
		deltamin = delta[t][survivorStates[t]];
		// s ... state aus dem man in den survivor state kommen h?tte k?nnen (2. M?glichkeit)

		if (previousMatrixColumn[t][survivorStates[t]] == -1)
		{
			// wenn keine Entscheidung getroffen wurde muss Bit flippen
			s = previousState(survivorStates[t], survivorBit[t][survivorStates[t]] * (-1) + 1);
		}
		else
		{
			// falls Entscheidung getroffen wird muss Spalte wechseln
			// 2 F?lle:	a) 0 <--> 2
			//			b) 1 <--> 2
			// ob 0 oder 1 steht im survivorBit
			int column;
			if (previousMatrixColumn[t][survivorStates[t]] == 2)
			{
				column = survivorBit[t][survivorStates[t]];
			}
			else
			{
				column = 2;
			}
			//s = previousState[previousMatrixColumn[t][survivorStates[t]]*(-1)+1][survivorStates[t]][survivorBit[t][survivorStates[t]]];
			s = previousState(survivorStates[t], column);
		}

		for (int i = t-1; i > 0; i--)
		{
			//sucht das kleinste delta auf dem survivor state
			if (delta[i][survivorStates[i]] < deltamin)
			{
				deltamin = delta[i][survivorStates[i]];
			}
			if (survivorBit[i][survivorStates[i]] != survivorBit[i][s])
			{
				delta[i][survivorStates[i]] = deltamin;
				#if DEBUG == 1
				printf("\nupdate k=%d i=%d, delta=%f",t,i,delta[i][survivorStates[i]]);
				#endif
			}

			if (previousMatrixColumn[i][s] == -1)
			{
				//wenn keine Entscheidung getroffen wurde muss Bit flippen
				s = previousState(s, survivorBit[i][s]);
			}
			else
			{
				//falls Entscheidung getroffen wird muss Spalte wechseln
				s = previousState(s, previousMatrixColumn[i][s]);
			}
		}

	}

	NumericVector softOutput(msgLen);

	for(int t = 1; t < msgLen+1; t++)
	{
		softOutput[t-1] = delta[t][survivorStates[t]] * (survivorBit[t][survivorStates[t]]*2-1) - La[t-1] - Lc * x_d[t-1];
	}

	#if DEBUG == 1
	printf("\n");
	for (int i = 0; i < msgLen; i++)
	{
		printf("%f\n",softOutput[i]);
	}
	#endif

	return softOutput;
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
	IntegerMatrix previousState,
	IntegerMatrix output,
	int output_index
)
{
	const int msgLen = x_noisy.size();

	NumericVector x_d_p(msgLen);  //noisy data permutated
	NumericVector Le1;    //decoder #1 extrinsic likelihood
	NumericVector Le1_p(msgLen);  //decoder #1 extrinsic likelihood permuted
	NumericVector Le2;    //decoder #2 extrinsic likelihood
	NumericVector Le2_ip(msgLen); //decoder #2 extrinsic likelihood inverse permuted


    //zero apriori information into very first iteration of BCJR
    for(int k = 0; k < msgLen; k++)
		Le2_ip[k] = 0;

    for(int i = 0; i < N_ITERATION; i++)
    {
		#if MAP_SOVA == 0
    	//modified_bcjr(Lc, Le2_ip, x_d, p1_d, Le1, 1);
		#elif MAP_SOVA == 1
		Le1 = c_sova(x_noisy, parity_noisy1, Le2_ip, 1, N, M, previousState, output, output_index);
		#endif

       //permute decoder#1 likelihoods to match decoder#2 order
    	for(int k = 0; k < msgLen; k++)
    	{
			Le1_p[k] = Le1[permutation[k]];
			x_d_p[k] = x_noisy[permutation[k]];
    	//	printf("%f ",x_d_p[k]=x_d[permutation[k]]);
    	//	printf("%f ",p2_d[k]);
    	//printf("\n");
    	}
		#if MAP_SOVA == 0
		//modified_bcjr(Lc, Le1_p,  x_d_p, p2_d, Le2, 0);
		#elif MAP_SOVA == 1
		Le2 = c_sova(x_noisy, parity_noisy2, Le1_p, 0, N, M, previousState, output, output_index);
		#endif

        //inverse permute decoder#2 likelihoods to match decoder#1 order
    	for(int k = 0; k < msgLen; k++)
		{
    		Le2_ip[permutation[k]] = Le2[k];
		}

        #if DEBUG > 0
		for(int k = 0; k < msgLen; k++)
		{
 			printf("\ni=%d Le1[%i]=%f\t",i+1, k, Le1[k]);
 			printf("Le2_ip[%i]=%f\t", k, Le2_ip[k]);
 			//printf("Le2[%i] = %f\t", k, Le2[k]);
 			//printf("Lc*x_d[%i] = %f", k, Lc*x_d[k]);
			printf("L[%i] = %f",k, Lc*x_noisy[k] + Le1[k] + Le2_ip[k]);
		}
		printf("\n");
		#endif
	}

	NumericVector softOutput(msgLen);
	IntegerVector hardOutput(msgLen);

    //calculate overall likelihoods and then slice'em
    for(int k = 0; k < msgLen; k++)
    {
		softOutput[k] = Lc * x_noisy[k] + Le1[k] + Le2_ip[k]; 		//soft decision
		hardOutput[k] = (softOutput[k] > 0.0) ? 1 : 0;         //hard decision
	}


    //print soft decisions
	for(int k = 0; k < msgLen ; k++)
		printf("L_h[%i] = %f\n", k, softOutput[k]);
	printf("\n");

    //print hard decisions
	printf("X_h = ");
    for(int k = 0; k < msgLen; k++)
    	printf("%i", hardOutput[k]);

    List result = List::create(Rcpp::Named("softOutput") = softOutput,
							   Rcpp::Named("hardOutput") = hardOutput);

	return result;
}
