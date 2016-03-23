#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAP_SOVA 1		//set which algorithm to use [0=>MAP,1=>SOVA]
#define DEBUG 1
//noise standard deviation
double sigma = 1.0;	//Lc=2/(sigma * sigma)=2

//
//       +-----------> Xk
//       |
//       |
//       |
//  Xk---+--[E1]-----> P1k
//       |
//      [P]
//       |
//       +--[E2]-----> P2k
//
//
void turbo_encode
(
	int N,					//length of information block
	int X[N],   			//block of N-2 information bits + 2 to_be_decided bits
	int P1[N],  			//encoder #1 parity bits
	int P2[N],   			//encoder #2 parity bits
	int permutation[N],		//interleaver array
	int M,					//amount of states
	int term_bits,			//amount of termination bits
	int next[M][2],
	int parity[M][2],
	int term[M][term_bits],
	int use_parity_index	//which output of parity should be used
)
{
	int	state;	//encoder state
	int	X_p[N];	//X_permuted = permuted bits 

	//encoder #1
	state = 0;		//encoder always starts in state 0
	for(int k = 0; k < N-term_bits; k++)
	{
		P1[k] = (parity[state][X[k]] >> use_parity_index) & 1;
		state = next[state][X[k]];
		//printf("s[%i] = %i\n", k, state);
	}

	int j = 0;
	for(int i = term_bits; i > 0; i--)
	{
		//terminate encoder #1 trellis to state 0
		X[N-i]  = term[state][j];  //databit to feed a 0 into delay line
		j++;
	}

	for(int i = term_bits; i > 0; i--)
	{
		P1[N-i] = (parity[state][X[N-i]] >> use_parity_index) & 1; //parity from state with databit included termination
		state   = next[state][X[N-i]];   //next state from current state
	}

	if(state != 0)
	{
		//should never get here
		printf("Error: Could not terminate encoder #1 trellis\n");
		exit(1);
	}

	//permute tx databits for encoder #2
	for(int k = 0; k < N; k++)
		X_p[k] = X[permutation[k]];

	//encoder #2
	state = 0; //encoder always starts in state 0
	for(int k = 0; k < N; k++)
	{
		P2[k] = (parity[state][X_p[k]] >> use_parity_index) & 1; //parity from state with databit=X_p[k]
		state = next[state][X_p[k]];   //next state from current state
	}

	//for(int k = 0; k < N; k++)
	//	printf("%i %i %i %i\n", X[k], P1[k], X_p[k], P2[k]);

}

//  SOVA algorithm (MAP decoder)
//
void sova
(
	int		N,				//so its resolving from state 0 , otherwise its resolving from max metric state
	double	Lc,				//Lc = 2/(sigma*sigma) = channel reliability
	double	La[N],			//apriori likelihood of each info bit
	double	x_d[N],			//noisy data
	double	p_d[N],			//noisy parity
	double	Le[N],			//extrinsic log likelihood
	int		IsTerminated,	//usualy first component decoder is terminated (this variable is 1) 
	int		M,					//amount of states
	int		term_bits,			//amount of termination bits
	int		previous[2][M][2],
	int		parity[M][2],
	int		use_parity_index	//which output of parity should be used
)
{
	double Metric[N+1][M];
	double Max[2];
	double Metric_temp;
	double delta[N+1][M];
	double deltamin;
	int survivor_bit[N+1][M];
	int survivor_states[N+1];
	int previous_matrix_decision[N+1][M];

	Metric[0][0]=0;
	for(int s = 1; s < M; s++)
	{
		Metric[0][s] = -1000;		
	}
	for(int i = 0; i < N+1; i++)
	{
		for(int j = 0; j < M; j++)
		{
			previous_matrix_decision[i][j] = -1;	
		}
	}

	for(int k = 1; k < N+1; k++)
	{
		for(int s = 0; s < M; s++)
		{
			for(int i = 0; i < 2; i++)
			{
				//when previous doesn't exists
				if(previous[0][s][i] != -1)
				{ 
					Max[i] = Metric[k-1][previous[0][s][i]]  +  (i * 2-1) * (Lc * x_d[k-1] + La[k-1]) + ((((parity[previous[0][s][i]][i] >> use_parity_index) & 1) > 0) * 2 -1) * Lc * p_d[k-1];	
					if(previous[1][s][i] != -1)
					{
						previous_matrix_decision[k][s] = 0;
						Metric_temp = Metric[k-1][previous[1][s][i]]  +  (i * 2-1) * (Lc * x_d[k-1] + La[k-1]) + ((((parity[previous[1][s][i]][i] >> use_parity_index) & 1) > 0) * 2 -1) * Lc * p_d[k-1];
						delta[k][s] = abs(Max[i] - Metric_temp) / 2;
						survivor_bit[k][s] = i;
						if(Metric_temp > Max[i])
						{
							Max[i] = Metric_temp;
							previous_matrix_decision[k][s] = 1;	
						}
						Metric[k][s] = Max[i];
					}
				} else
				{
					Max[i] = -1000;
				}
			}
						
			//wird nur ausgeführt, wenn keine Entscheidung getroffen werden musste
			if(previous_matrix_decision[k][s] == -1)
			{
				if (Max[0] > Max[1])
				{
					survivor_bit[k][s] = 0;
					Metric[k][s] = Max[0];
				} else
				{
					survivor_bit[k][s] = 1;
					Metric[k][s] = Max[1];
				}
				
				delta[k][s] = abs(Max[0]- Max[1]) / 2;
			}
			
			
			/* Belegung previous_matrix_decision
			*  -1, wenn keine Entscheidung getroffen wurde (Entscheidung zwischen beiden Matrizen)
			*  0,1 die jeweilige Entscheidung, ob Matrix 0 oder 1 verwendet wird (1.Index der previous-Matrix)
			*/
			
			#if DEBUG == 2
			printf("\ndelta[%d][%d]=%f\t M0=%.2f \tM1=%.2f \tmax=%f",k,s,delta[k][s], Max[0],Max[1], Metric[k][s]);
			#endif
		}
//		printf("\n");
	}

	//s = previous[0][survivor[N][0]]; //predhodni iz stanja 0 po bitu (preziveli u N tom koraku ka stanju 0)
	if (IsTerminated > 0) 
	{
		survivor_states[N] = 0;	
	}
	else
	{
		double max_metric = Metric[N][0];
		survivor_states[N] = 0;
		for(int s = 1; s < M; s++)
		{	
			if (max_metric < Metric[N][s])
			{
				max_metric = Metric[N][s];
				survivor_states[N] = s;
			}
		}
	}

	for (int k = N-1; k >= 0; k--)
	{
		if(previous_matrix_decision[k+1][survivor_states[k+1]] == -1)
		{
			survivor_states[k] = previous[0][survivor_states[k+1]][survivor_bit[k+1][survivor_states[k+1]]];
		} else 
		{
			survivor_states[k] = previous[previous_matrix_decision[k+1][survivor_states[k+1]]][survivor_states[k+1]][survivor_bit[k+1][survivor_states[k+1]]];
		}
		
		if(survivor_states[k] == -1){
			printf("survivor_state -1!!!\n");
			exit(0);
		}

		//printf("\nk=%d s=%d, delta=%f, b=%d", k , survivor_states[k],delta[k][survivor_states[k]],survivor_bit[k][survivor_states[k]] );
	}

	int s;
	for (int k = term_bits+1; k < N+1; k++)
	{
		deltamin = delta[k][survivor_states[k]];
		// s ... state aus dem man in den survivor state kommen hätte können (2. Möglichkeit)
		if(previous_matrix_decision[k][survivor_states[k]] == -1)
		{
			//wenn keine Entscheidung getroffen wurde muss Bit flippen
			s = previous[0][survivor_states[k]][survivor_bit[k][survivor_states[k]] * (-1) + 1];
		} else
		{
			//falls Entscheidung getroffen wird muss State flippen
			s = previous[previous_matrix_decision[k][survivor_states[k]]*(-1)+1][survivor_states[k]][survivor_bit[k][survivor_states[k]]];
		}

		for (int i = k-1; i > 0; i--)
		{
			//sucht das kleinste delta auf dem survivor state
			if (delta[i][survivor_states[i]] < deltamin)
			{
				deltamin=delta[i][survivor_states[i]];				
			}
			
			//wenn kein voriger Zustand in previous existiert
			if(s == -1)
			{
				printf("Bei delta-Berechnung, survivor_state -1!!!\n");
				exit(0);
			}

			if (survivor_bit[i][survivor_states[i]] != survivor_bit[i][s])
			{
			//	if (delta[i][survivor_states[i]] != deltamin)
				//	{
					delta[i][survivor_states[i]] = deltamin;
					#if DEBUG==2
					printf("\nupdate k=%d i=%d, delta=%f ", k,i, delta[i][survivor_states[i]]);
					#endif
			//	}
			}
			
			if(previous_matrix_decision[i][s] == -1)
			{
				//wenn keine Entscheidung getroffen wurde muss Bit flippen
				s = previous[0][s][survivor_bit[i][s]];
			} else
			{
				//falls Entscheidung getroffen wird muss State flippen
				s = previous[previous_matrix_decision[i][s]][s][survivor_bit[i][s]];
			}
		}

	}

	for(int k = 1; k < N+1; k++)
	{
//		printf("\ndelta[%d]=%f",k,delta[k][survivor_states[k]]);
		Le[k-1]= delta[k][survivor_states[k]] * (survivor_bit[k][survivor_states[k]]*2-1) - La[k-1] - Lc * x_d[k-1];
//		printf("\nLe[%d]=%f ", k-1,Le[k-1]);
	}
}

void turbo_decode(
	int		N,					//length of information block
	double	sigma,   			//channel noise standard deviation
	double	x_d[N],  			//x_dash  = noisy data symbol   
	double	p1_d[N], 			//p1_dash = noisy parity#1 symbol
	double	p2_d[N], 			//p2_dash = noisy parity#2 symbol
	double	L_h[N],  			//L_hat = likelihood of databit given entire observation
	int		X_h[N],  			//X_hat = sliced MAP decisions
	int		permutation[N],		//interleaver array
	int		M,					//amount of states
	int		term_bits,			//amount of termination bits
	int		previous[2][M][2],
	int		next[M][2],
	int		parity[M][2],
	int		term[M][term_bits],
	int		use_parity_index,	//which output of parity should be used
	int		N_ITERATION			//amount of iterations by decoding process
)
{
	double x_d_p[N];  //noisy data permutated
	double Le1[N];    //decoder #1 extrinsic likelihood
	double Le1_p[N];  //decoder #1 extrinsic likelihood permuted
	double Le2[N];    //decoder #2 extrinsic likelihood
	double Le2_ip[N]; //decoder #2 extrinsic likelihood inverse permuted
    double Lc;        //channel reliability value

    Lc = 2.0 / (sigma*sigma); //requires sigma to be non-trivial

    //zero apriori information into very first iteration of BCJR
    for(int k = 0; k < N; k++)
		Le2_ip[k] = 0;

    for(int i = 0; i < N_ITERATION; i++)
    {
		#if MAP_SOVA == 0 
    	modified_bcjr(Lc, Le2_ip, x_d, p1_d, Le1, 1);
		#elif MAP_SOVA == 1
		sova(N, Lc, Le2_ip, x_d, p1_d, Le1, 1, M, term_bits, previous, parity, use_parity_index);
		#endif

       //permute decoder#1 likelihoods to match decoder#2 order
    	for(int k = 0; k < N; k++)
    	{	
			Le1_p[k] = Le1[permutation[k]];
			x_d_p[k] = x_d[permutation[k]];
    //	printf("%f ",x_d_p[k]=x_d[permutation[k]]);
    //	printf("%f ",p2_d[k]);
    	//printf("\n");
    	}
		#if MAP_SOVA == 0 
		modified_bcjr(Lc, Le1_p,  x_d_p, p2_d, Le2, 0);
		#elif MAP_SOVA == 1
		sova(N, Lc, Le1_p,  x_d_p, p2_d, Le2, 0, M, term_bits, previous, parity, use_parity_index);
		#endif

        //inverse permute decoder#2 likelihoods to match decoder#1 order
    	for(int k = 0; k < N; k++)
		{
    		Le2_ip[permutation[k]] = Le2[k];	
		}

        #if DEBUG > 0
		for(int k = 0; k < N; k++)
		{
 			printf("\ni=%d Le1[%i]=%f\t",i+1, k, Le1[k]);
 			printf("Le2_ip[%i]=%f\t", k, Le2_ip[k]);
 			//printf("Le2[%i] = %f\t", k, Le2[k]);
 			//printf("Lc*x_d[%i] = %f", k, Lc*x_d[k]);
			printf("L[%i] = %f",k, Lc*x_d[k] + Le1[k] + Le2_ip[k]);
		}
		printf("\n");
		#endif
	}

    //calculate overall likelihoods and then slice'em
    for(int k = 0; k < N; k++)
    {
		L_h[k] = Lc*x_d[k] + Le1[k] + Le2_ip[k]; //soft decision
		X_h[k] = (L_h[k] > 0.0) ? 1 : 0;         //hard decision
	}
}

void wrapper_encode
(
	int *input,
	int *output,
	int *permutation,
	int *input_length,
	int *output_length,
	int *amount_register,
	int *previous_in,
	int *next_in,
	int *parity_in,
	int *term_in,
	int *use_parity_index
)
{
	int    parity_out1[*input_length];     //encoder #1 parity bits
	int    parity_out2[*input_length];     //encoder #2 parity bits
		
	int amount_states = pow(2, *amount_register);
	int previous[2][amount_states][2];   //previous[m][i] = previous state (this comes to state m with databit = i)
	int next[amount_states][2];     //next[m][i] = next state (this comes from state m with databit = i)
	int parity[amount_states][2]; //parity bit associated with transition from state m
	int term[amount_states][*amount_register];   //term[m] = pair of data bits required to terminate trellis
	
	for(int i = 0; i < amount_states; i++)
	{
		for(int j = 0; j < 2; j++)
		{ 
			next[i][j] = next_in[i*2+j]; 
			parity[i][j] = parity_in[i*2+j];  
		}
	}
	for(int i = 0; i < amount_states; i++)
	{
		for(int j = 0; j < *amount_register; j++)
		{
			term[i][j] = term_in[i*(*amount_register)+j];  
		}
	}
	
	//i = Matrixauswahl, j = State, i = Bit
	for(int i = 0; i < 2; i++)
	{
		for(int j = 0; j < amount_states; j++)
		{
			for(int k = 0; k < 2; k++)
			{
				previous[i][j][k] = previous_in[amount_states*2*i+j*2+k];
			}
		}	
	}

		
	 /********************************
     *           ENCODER            *
     ********************************/

    turbo_encode(*input_length, input, parity_out1, parity_out2, permutation, amount_states, *amount_register, next, parity, term, *use_parity_index);

	//map databit to symbols and concat them to the output array
	for(int i = 0; i < *output_length; i++)
	{
		if(i < *input_length)
		{
			output[i] = input[i] ? +1 : -1;			
		}

		if(i >= *input_length && i < 2*(*input_length))
		{
			output[i] = parity_out1[i- (*input_length)] ? +1 : -1;			
		}

		if(i >= 2*(*input_length))
		{
			output[i] = parity_out2[i- 2*(*input_length)] ? +1 : -1;	
		}
	}
}

void wrapper_decode
(
	double *input,
	double *output_soft,
	int *output_hard,
	int *permutation,
	int *input_length,
	int *output_length,
	int *iterations,
	int *amount_register,
	int *previous_in,
	int *next_in,
	int *parity_in,
	int *term_in,
	int *use_parity_index
)
{
	double x_noisy[*output_length];    		//x_dash  = noisy data symbol
	double parity_noisy1[*output_length];   //p1_dash = noisy parity#1 symbol
	double parity_noisy2[*output_length];   //p2_dash = noisy parity#2 symbol
	
	int amount_states = pow(2, *amount_register);
	int previous[2][amount_states][2];   //previous[m][i] = previous state (this comes to state m with databit = i)
	int next[amount_states][2];     //next[m][i] = next state (this comes from state m with databit = i)
	int parity[amount_states][2]; //parity bit associated with transition from state m
	int term[amount_states][*amount_register];   //term[m] = pair of data bits required to terminate trellis
	
	for(int i = 0; i < amount_states; i++)
	{
		for(int j = 0; j < 2; j++)
		{
			next[i][j] = next_in[i*2+j]; 
			parity[i][j] = parity_in[i*2+j]; 
		}
	}
	
	for(int i = 0; i < amount_states; i++)
	{
		for(int j = 0; j < *amount_register; j++)
		{
			term[i][j] = term_in[i*(*amount_register)+j];  
		}
	}
	
	//i = Matrixauswahl, j = State, i = Bit
	for(int i = 0; i < 2; i++)
	{
		for(int j = 0; j < amount_states; j++)
		{
			for(int k = 0; k < 2; k++)
			{
				previous[i][j][k] = previous_in[amount_states*2*i+j*2+k];
			}
		}	
	}
	
	
	for(int i = 0; i < *input_length; i++)
	{
		if(i < *output_length)
		{
			x_noisy[i] = input[i];
		}
		if(i >= *output_length && i < 2*(*output_length))
		{
			parity_noisy1[i- (*output_length)] = input[i];	
		}
		if(i >= 2*(*output_length))
		{
			parity_noisy2[i- 2*(*output_length)] = input[i];	
		}
	}
		
	 /********************************
     *           DECODER            *
     ********************************/

    turbo_decode(	*output_length, 
					1,
					x_noisy,
					parity_noisy1,
					parity_noisy2,
					output_soft,
					output_hard,
					permutation,
					amount_states,
					*amount_register,
					previous,
					next,
					parity,
					term,
					*use_parity_index,
					*iterations	);
    
    //print soft decisions
	for(int k = 0; k < *output_length ; k++)
		printf("L_h[%i] = %f\n", k, output_soft[k]);
	printf("\n");

    //print hard decisions
	printf("X_h = ");
    for(int k = 0; k < *output_length; k++)
    	printf("%i", output_hard[k]);
}

void main(){
	
	double x_d[] = {1,-1,1,1,1,-1,-1,-1};
	double p1_d[] = {1,1,-1,-1,1,-1,1,-1};
	double Le1[8];
	double Le2_ip[8];
	for(int k = 0; k < 8; k++)
		Le2_ip[k] = 0;
		
	int previous[2][4][2];
	int parity[4][2];
	
	previous[0][0][0] = 0;
	previous[0][0][1] = -1; 
	previous[0][1][0] = 2;
	previous[0][1][1] = -1;
	previous[0][2][0] = -1;
	previous[0][2][1] = 0;
	previous[0][3][0] = -1;
	previous[0][3][1] = 2;
	previous[1][0][0] = 1;
	previous[1][0][1] = -1;
	previous[1][1][0] = 3;
	previous[1][1][1] = -1;
	previous[1][2][0] = -1;
	previous[1][2][1] = 1;
	previous[1][3][0] = -1;
	previous[1][3][1] = 3;
	
	parity[0][0] = 0;
	parity[0][1] = 3;
	parity[1][0] = 3;
	parity[1][1] = 0;
	parity[2][0] = 1;
	parity[2][1] = 3;
	parity[3][0] = 1;
	parity[3][1] = 3;
		
	sova(8, 2, Le2_ip, x_d, p1_d, Le1, 1, 4, 2, previous, parity, 0);
    
    //print soft decisions
	for(int k = 0; k < 8 ; k++)
		printf("\nL_h[%i] = %f", k, Le1[k]);
	printf("\n");
}
