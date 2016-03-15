#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAP_SOVA 1		//set which algorithm to use [0=>MAP,1=>SOVA]
    //noise standard deviation
double sigma = 1.0;	//Lc=2/(sigma * sigma)=2

// tables for trellis are global
//
#define M  4      //no. of trellis states

int previous[M][2];   //previous[m][i] = previous state (this comes to state m with databit = i)
int next[M][2];     //next[m][i] = next state (this comes from state m with databit = i)
int parity[M][2]; //parity bit associated with transition from state m
int term[M][2];   //term[m] = pair of data bits required to terminate trellis

//
//      +--------------------------> Xk
//      |  fb
//      |  +---------(+)-------+
//      |  |          |        |
//  Xk--+-(+)-+->[D]----->[D]--+
//            |                |
//            +--------------(+)---> Pk
//
//
void gen_tab(void)
{
	int m, i, b0, b1, fb, state;

    //generate tables for 4 state RSC encoder
	for(m = 0; m < M; m++) //for each starting state
	{
		for(i = 0; i < 2; i++) //for each possible databit
		{
			b0 = (m >> 0) & 1; //bit 0 of state
			b1 = (m >> 1) & 1; //bit 1 of state

			//parity from state m with databit i
			parity[m][i] = b0 ^ i;
			
			//next[m][i] = next state from state m with databit i
			next[m][i]   = b0 *2 + (i ^ b0 ^ b1);
		//	printf("next[%d][%d]=%d\n",m,i,next[m][i]);
		}
	//	printf("parity[%d]= %d %d\n",m,parity[m][0],parity[m][1]);
	}

    //from[m][i] = previous state to state m with databit i
    for(m = 0; m < M; m++)
    	for(i = 0; i < 2; i++)
			previous[next[m][i]][i] = m;

	//  Generate table of data bit pairs which terminate
	//  the trellis for a given state m
	//
	//  We simply set Xk equal to the feedback value to force
	//  the delay line to fill up with zeros.
	//
	for(m = 0; m < M; m++) //for each state
	{
		state = m;
		b0 = (state >> 0) & 1; //bit 0 of state
		b1 = (state >> 1) & 1; //bit 1 of state
		fb = b0 ^ b1;          //feedback bit
		term[m][0] = fb;       //will set X[N-2] = fb

		state = next[m][fb];   //advance from state m with databit=fb
		b0 = (state >> 0) & 1; //bit 0 of state
		b1 = (state >> 1) & 1; //bit 1 of state
		fb = b0 ^ b1;          //feedback bit
		term[m][1] = fb;       //will set X[N-1] = fb
	}
}

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
	int N,		//length of information block
	int X[N],   //block of N-2 information bits + 2 to_be_decided bits
	int P1[N],  //encoder #1 parity bits
	int P2[N],   //encoder #2 parity bits
	int permutation[N]	//interleaver array
)
{
	int    k;      //trellis stage
	int    state;  //encoder state
	int    X_p[N]; //X_permuted = permuted bits

	//encoder #1
	state = 0; //encoder always starts in state 0
	for(k = 0; k < N-2; k++)
	{
		P1[k] = parity[state][X[k]];
		state = next[state][X[k]];
		//printf("s[%i] = %i\n", k, state);
	}

	//terminate encoder #1 trellis to state 0
	X[N-2]  = term[state][0];  //databit to feed a 0 into delay line
	X[N-1]  = term[state][1];  //databit to feed another 0 into delay line

	P1[N-2] = parity[state][X[N-2]]; //parity from state with databitX[N-2]
	state   = next[state][X[N-2]];   //next state from current state
    P1[N-1] = parity[state][X[N-1]]; //parity from state with databit=X[N-1]
	state   = next[state][X[N-1]];   //next state from current state

	if(state != 0)
	{
		//should never get here
		printf("Error: Could not terminate encoder #1 trellis\n");
		exit(1);
	}

	//permute tx databits for encoder #2
	for(k = 0; k < N; k++)
		X_p[k] = X[permutation[k]];

	//encoder #2
	state = 0; //encoder always starts in state 0
	for(k = 0; k < N; k++)
	{
		P2[k] = parity[state][X_p[k]]; //parity from state with databit=X_p[k]
		state = next[state][X_p[k]];   //next state from current state
	}

	//for(k = 0; k < N; k++)
//		printf("%i %i %i %i\n", X[k], P1[k], X_p[k], P2[k]);

}

//  SOVA algorithm (MAP decoder)
//
void sova
(
	int    N,				//so its resolving from state 0 , otherwise its resolving from max metric state
	double Lc,           //Lc = 2/(sigma*sigma) = channel reliability
	double La[N],        //apriori likelihood of each info bit
	double x_d[N],       //noisy data
	double p_d[N],       //noisy parity
	double Le[N],         //extrinsic log likelihood
	int    IsTerminated	//usualy first component decoder is terminated (this variable is 1) 
)
{
	int i,s,k;
	double Metric[N+1][M];
	double Max[2];
	double delta[N+1][M];
	double deltamin;
	int survivor_bit[N+1][M]; //bit koji je preziveo
	int survivor_states[N+1];

//	for(k=0;k<N;k++)
//		printf("\n%f\t%f\t%f",La[k],x_d[k],p_d[k]);

	Metric[0][0]=0;
	for(s=1;s<M;s++)
		Metric[0][s]=-1000;

	for(k=1;k<N+1;k++)
	{
		for(s=0;s<M;s++)
		{
			for(i=0;i<2;i++)
			{
				Max[i]=Metric[k-1][previous[s][i]]  +  (i*2-1) * (Lc * x_d[k-1] + La[k-1]) + ((parity[previous[s][i]][i]>0) * 2 -1)*Lc * p_d[k-1];
			}
			if (Max[0]>Max[1])
			{
				survivor_bit[k][s]=0;
				Metric[k][s]=Max[0];
				delta[k][s]=(Max[0]-Max[1])/2;
			}
			else
			{
				survivor_bit[k][s]=1;
				Metric[k][s]=Max[1];
				delta[k][s]=(Max[1]-Max[0])/2;
			}
			#if DEBUG == 2
			printf("\ndelta[%d][%d]=%f\t M0=%.2f \tM1=%.2f \tmax=%f",k,s,delta[k][s], Max[0],Max[1], Metric[k][s]);
			#endif
		}
//		printf("\n");
	}

	//sad se vracamo po pobednickoj putanji iz stanja 0, i zapamtpimo pobednicku putanju
	//zatim idemo od pocetka po pobednickoj putanji i pamtimo pozicije gde se razlikuju pobednicka i izgubljena putanja
	//jer samo te delte updejtujemo

	//s = previous[0][survivor[N][0]]; //predhodni iz stanja 0 po bitu (preziveli u N tom koraku ka stanju 0)
	if (IsTerminated>0) 
	{
		survivor_states[N]=0;	
	}
	else
	{
		double max_metric = Metric[N][0];
		survivor_states[N]=0;
		for(s=1;s<M;s++)
		{	
			if (max_metric < Metric[N][s])
			{
				max_metric = Metric[N][s];
				survivor_states[N]=s;
			}
		}
	}

	for (k=N-1;k>=0;k--)
	{
		survivor_states[k] = previous[survivor_states[k+1]][survivor_bit[k+1][survivor_states[k+1]]];
//		printf("\nk=%d s=%d, delta=%f, b=%d", k , survivor_states[k],delta[k][survivor_states[k]],survivor_bit[k][survivor_states[k]] );
	}
	for (k=3;k<N+1;k++)
	{
		deltamin = delta[k][survivor_states[k]];
		// s ... state aus dem man in den survivor state kommen hätte können (2. Möglichkeit)
		s = previous[ survivor_states[k] ][ survivor_bit[k][survivor_states[k] ] * (-1) + 1]; //competing putanja
		for (i=k-1;i>0;i--)
		{
			//sucht das kleinste delta auf dem survivor state
			if (delta[i][survivor_states[i]] < deltamin)
								deltamin=delta[i][survivor_states[i]];
			if (survivor_bit[i][survivor_states[i]]!=survivor_bit[i][s])
			{
			//	if (delta[i][survivor_states[i]] != deltamin)
				//	{
					delta[i][survivor_states[i]] = deltamin;
					#if DEBUG==2
					printf("\nupdate k=%d i=%d, delta=%f ", k,i, delta[i][survivor_states[i]]);
					#endif
			//	}
			}
			s = previous[s][survivor_bit[i][s]];
		}

	}

	for(k=1;k<N+1;k++)
	{
//		printf("\ndelta[%d]=%f",k,delta[k][survivor_states[k]]);
		Le[k-1]= delta[k][survivor_states[k]] * (survivor_bit[k][survivor_states[k]]*2-1) - La[k-1] - Lc * x_d[k-1];
//		printf("\nLe[%d]=%f ", k-1,Le[k-1]);
	}
}

void turbo_decode(
	int	   N,				//length of information block
	double sigma,   //channel noise standard deviation
	double x_d[N],  //x_dash  = noisy data symbol   
	double p1_d[N], //p1_dash = noisy parity#1 symbol
	double p2_d[N], //p2_dash = noisy parity#2 symbol
	double L_h[N],  //L_hat = likelihood of databit given entire observation
	int    X_h[N],  //X_hat = sliced MAP decisions
	int    permutation[N],	//interleaver array
	int	   N_ITERATION		//amount of iterations by decoding process
)
{
	int i, k;

	double x_d_p[N];  //
	double Le1[N];    //decoder #1 extrinsic likelihood
	double Le1_p[N];  //decoder #1 extrinsic likelihood permuted
	double Le2[N];    //decoder #2 extrinsic likelihood
	double Le2_ip[N]; //decoder #2 extrinsic likelihood inverse permuted
    double Lc;        //channel reliability value

    Lc = 2.0 / (sigma*sigma); //requires sigma to be non-trivial

    //zero apriori information into very first iteration of BCJR
    for(k = 0; k < N; k++)
		Le2_ip[k] = 0;

    for(i = 0; i < N_ITERATION; i++)
    {
		#if MAP_SOVA == 0 
    	modified_bcjr(Lc, Le2_ip, x_d, p1_d, Le1, 1);
		#elif MAP_SOVA == 1
		sova(N, Lc, Le2_ip, x_d, p1_d, Le1, 1);
		#endif

       //permute decoder#1 likelihoods to match decoder#2 order
    	for(k = 0; k < N; k++)
    	{	
			Le1_p[k] = Le1[permutation[k]];
			x_d_p[k]=x_d[permutation[k]];
    //	printf("%f ",x_d_p[k]=x_d[permutation[k]]);
    //	printf("%f ",p2_d[k]);
    	//printf("\n");
    	}
		#if MAP_SOVA == 0 
		modified_bcjr(Lc, Le1_p,  x_d_p, p2_d, Le2, 0);
		#elif MAP_SOVA == 1
		sova(N, Lc, Le1_p,  x_d_p, p2_d, Le2, 0);
		#endif

        //inverse permute decoder#2 likelihoods to match decoder#1 order
    	for(k = 0; k < N; k++)
    		Le2_ip[permutation[k]] = Le2[k];

        #if DEBUG > 0
		for(k = 0; k < N; k++)
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
    for(k = 0; k < N; k++)
    {
		L_h[k] = Lc*x_d[k] + Le1[k] + Le2_ip[k]; //soft decision
		X_h[k] = (L_h[k] > 0.0) ? 1 : 0;         //hard decision
	}
}

//// A utility function to swap to integers
//void swap (int *a, int *b)
//{
//    int temp = *a;
//    *a = *b;
//    *b = temp;
//}
//
//void create_permutation_array(int *permutation, int N){
//
//	for(int i = 0; i < N; i++){
//		permutation[i] = i;
//	}
//	for (int i = N-1; i > 0; i--){
//        // Pick a random index from 0 to i
//        int j = rand() % (i+1);
// 
//        // Swap arr[i] with the element at random index
//        swap(&permutation[i], &permutation[j]);
//    }
//}

void wrapper_encode(int *input, int *output, int *permutation, int *input_length, int *output_length){
	int    P1[*input_length];     //encoder #1 parity bits
	int    P2[*input_length];     //encoder #2 parity bits
	
	srand(1);    //init random number generator
	
	gen_tab();   //generate trellis tables
	
	 /********************************
     *           ENCODER            *
     ********************************/

    turbo_encode(*input_length, input, P1, P2, permutation);

	for(int i = 0; i < *output_length; i++){
		if(i < *input_length)
			output[i] = input[i];
		if(i >= *input_length && i < 2*(*input_length))
			output[i] = P1[i- (*input_length)];
		if(i >= 2*(*input_length))
			output[i] = P2[i- 2*(*input_length)];
	}
}

void wrapper_decode(double *input, double *output_soft, int *output_hard,  int *permutation, int *input_length, int *output_length, int *iterations){
	double x_d[*output_length];    //x_dash  = noisy data symbol
	double p1_d[*output_length];   //p1_dash = noisy parity#1 symbol
	double p2_d[*output_length];   //p2_dash = noisy parity#2 symbol
	
	for(int i = 0; i < *input_length; i++){
		if(i < *output_length){
			x_d[i] = input[i];
		}
		if(i >= *output_length && i < 2*(*output_length)){
			p1_d[i- (*output_length)] = input[i];	
		}
		if(i >= 2*(*output_length)){
			p2_d[i- 2*(*output_length)] = input[i];	
		}
	}
	
	srand(1);    //init random number generator
	
	 /********************************
     *           DECODER            *
     ********************************/

    turbo_decode(*output_length, 1, x_d, p1_d, p2_d, output_soft, output_hard, permutation, *iterations);
    
    //print soft decisions
	for(int k = 0; k < *output_length ; k++)
		printf("L_h[%i] = %f\n", k, output_soft[k]);
	printf("\n");

    //print hard decisions
	printf("X_h = ");
    for(int k = 0; k < *output_length; k++)
    	printf("%i", output_hard[k]);
}
