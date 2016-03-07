// Author: Dusan Orlovic
// Previous Author:  VA
// Email:   contact@vashe.org
// Version: 0.1
// File:    turbo_example.c
//
// Rate 1/3 turbo coder simple example (for educational purposes)
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// References:
//
//   1.S.A. Barbulescu, "Iterative Decoding of Turbo Codes and Other
//     Concatenated Codes"
//     http://www.itr.unisa.edu.au/~steven/thesis/sab.ps.gz
//
//   2.W.E. Ryan,"A Turbo Code Tutorial" New Mexico State University
//     http://vada.skku.ac.kr/ClassInfo/digital-com2000/slides/turbo2c.pdf
//
//   3.M.C. Reed, S.S. Pietrobon, "Turbo-code termination schemes and a
//     novel alternative for short frames" IEEE Int. Symp. on Personal,
//     Indoor and Mobile Radio Commun., Taipei, Taiwan, pp. 354-358,
//     Oct. 1996.
//     http://www.itr.unisa.edu.au/~steven/turbo/PIMRC96.ps.gz
//
// Revisions: 0.1 fixed typo - La missing from modified_bcjr args
//
//
// Author:  dusan orlovic
// Email:   orlovic@uns.ac.rs
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define DEBUG     0   //set level of debug visibility [0=>off,1=>min,2=>max]
#define NOISEOFF  2   //set to suppress noise in channel [0=>off, 1=>read from file, 2=>awgn]
#define noise_file "awgn.txt" //filename where noise date is stored (each value in separate line)
#define N_ITERATION 5 //no. of turbo decoder iterations
#define MAP_SOVA 1		//set which algorithm to use [0=>MAP,1=>SOVA]
    //noise standard deviation
double sigma = 1.0;	//Lc=2/(sigma * sigma)=2

//  Block length in our example (6 databits + 2 termination bits)
//  A practical turbo code must use a much larger block length to
//  approach the Shannon limit.
#define N    8

//  Transmit bits for our simple example
    //  The two last bits will be inserted later by encoder#1
    //  to park it at state 0
    //
int X[N] = {1,1,0,0,1,0,-1,-1}; 

//  Pseudo random permutation (of length = the block length).
//
int permutation[N] = {7,5,1,2,4,3,6,0};


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

//  Normally distributed number generator (ubiquitous Box-Muller method)
//
double normal(void)
{
	double x, y, rr, randn;
	do{
        x  = (double) 2*rand()/RAND_MAX - 1.0; //uniform in range [-1,1]
        y  = (double) 2*rand()/RAND_MAX - 1.0; //uniform in range [-1,1]
        rr = x*x + y*y;
    } while( rr >= 1 );
    randn = x*sqrt((-2.0*log(rr))/rr);
  return(randn);
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
	int X[N],   //block of N-2 information bits + 2 to_be_decided bits
	int P1[N],  //encoder #1 parity bits
	int P2[N]   //encoder #2 parity bits
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

//  modified BCJR algorithm (MAP decoder)
//
void modified_bcjr
(
	double Lc,           //Lc = 2/(sigma*sigma) = channel reliability
	double La[N],        //apriori likelihood of each info bit
	double x_d[N],       //noisy data
	double p_d[N],       //noisy parity
	double Le[N],         //extrinsic log likelihood
	int    is_term      //indicates if trellis terminated
)
{
	int    k, m, i;
	double xk_h, pk_h;      //databit & parity associated with a branch
	double gamma[N][M][2];  //gammas for total likelihoods
	double gammae[N][M][2]; //gammas for extrinsic likelihoods
	double pr0, pr1;        //extrinsic likelihood = pr1/pr0
	double alpha[N+1][M];   //probability of entering branch via state m
	double beta[N+1][M];    //probability of exiting branch via state m
	double total;           //used for normalising alpha's and beta's

    //calculate branch gamma's
    for(k = 0; k < N; k++) //for each trellis stage
    {
		for(m = 0; m < M; m++) //for each state
		{
			for(i = 0; i < 2; i++) //for each databit
			{
				//data associated with branch
				xk_h = i ? +1 : -1;            //map databit to PAM symbol

				//parity associated with branch
				pk_h = parity[m][i] ? +1 : -1; //map parity bit to PAM symbol

                //used later to calculate alpha's and beta's
				gamma[k][m][i] = exp(0.5*(La[k] * xk_h +
				                          Lc * x_d[k] * xk_h +
				                          Lc * p_d[k] * pk_h));

                //used later to calculate extrinsic likelihood
				gammae[k][m][i] = exp(0.5*(Lc * p_d[k] * pk_h));
			}
		}
	}

	//  Calculate state alpha's
	//
    //  As the likelihood ratio for each stage k has a linear combination
    //  of alpha terms in both the numerator and the denominator, we can
    //  scale all the alpha's by any convenient scaling constant.
    //
    //  To help avoid underflow/overflow we normalise the alphas at each
    //  stage so that sum across all states is unity.
    //

    //  The encoders always start in state 0
	alpha[0][0] = 1;
	for(m = 1; m < M; m++)
		alpha[0][m] = 0;

	for(k = 1; k <= N; k++)
	{
		total = 0;

	    for(m = 0; m < M; m++)
	    {
			alpha[k][m] = alpha[k-1][previous[m][0]] * gamma[k-1][previous[m][0]][0] +
			              alpha[k-1][previous[m][1]] * gamma[k-1][previous[m][1]][1];

			total += alpha[k][m];
		}

		//normalise
		for(m = 0; m < M; m++)
			alpha[k][m] /= total;
	}

	//  Calculate state beta's
	//
    //  As the likelihood ratio for each stage k has a linear combination
    //  of beta terms in both the numerator and the denominator, we can
    //  scale all the beta's by any convenient scaling constant.
    //
    //  To help avoid underflow/overflow we normalise the betas at each
    //  stage so that sum across all states is unity.
    //

	if(is_term)                 //if trellis terminated
	{
		//we know for sure the final state is 0
	    beta[N][0] = 1;
	    for(m = 1; m < M; m++)
	    	beta[N][m] = 0;
	}
	else                       //else trellis not terminated
	{
		//we haven't a clue which is final state
		//so the best we can do is say they're all equally likely
	    for(m = 0; m < M; m++)
	    	beta[N][m] = 1.0 / (double) M;
	}

    //iterate backwards through trellis
	for(k = N-1; k >= 0; k--)
	{
		total = 0;
		for(m = 0; m < 4; m++)
		{
			beta[k][m] = beta[k+1][next[m][0]] * gamma[k][m][0] +
				         beta[k+1][next[m][1]] * gamma[k][m][1];


			total += beta[k][m];
		}

        //normalise
		for(m = 0; m < 4; m++)
			beta[k][m] /= total;
	}

    //  Calculate extrinsic likelihood
    //
	//  This is the information "gleaned" from the parity check
	//  Note the Ck's in equation 20 in [2] are different in the
	//  numerator and denominator. This is why the channel and
	//  apriori likelihoods can be brought out.
	//
	for(k = 0; k < N; k++)
	{
		pr0 = pr1 = 0;
		for(m = 0; m < 4; m++)
		{
			//we use gammae rather than gamma as we want the
			//extrinsic component of the overall likelihood
			pr1 += (alpha[k][m] * gammae[k][m][1] * beta[k+1][next[m][1]]);
			pr0 += (alpha[k][m] * gammae[k][m][0] * beta[k+1][next[m][0]]);
		}
		Le[k] = log(pr1 / pr0); //extrinsic likelihood
	}

    #if DEBUG > 1
    for(k = 0; k < N; k++)
    {
		for(m = 0; m < M; m++)
		{
			for(i = 0; i < 2; i++)
			{
				printf("gamma[%i][%i][%i]  = %f\t", k, m, i, gamma[k][m][i]);
				printf("gammae[%i][%i][%i] = %f\n", k, m, i, gammae[k][m][i]);
			}
		}
		printf("\n");
	}

	for(k = 0; k <= N; k++)
	{
	    for(m = 0; m < M; m++)
			printf("alpha[%i][%i] = %f\n", k, m, alpha[k][m]);
		printf("\n");
	}
	for(k = 0; k <= N; k++)
	{
	    for(m = 0; m < M; m++)
			printf("beta[%i][%i] = %f\n", k, m, beta[k][m]);
		printf("\n");
	}
    #endif

}

//  SOVA algorithm (MAP decoder)
//
void sova
(
	double Lc,           //Lc = 2/(sigma*sigma) = channel reliability
	double La[N],        //apriori likelihood of each info bit
	double x_d[N],       //noisy data
	double p_d[N],       //noisy parity
	double Le[N],         //extrinsic log likelihood
	int IsTerminated	//usualy first component decoder is terminated (this variable is 1) 
						//so its resolving from state 0 , otherwise its resolving from max metric state
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
		s = previous[ survivor_states[k] ][ survivor_bit[k][survivor_states[k] ] * (-1) + 1]; //competing putanja
		for (i=k-1;i>0;i--)
		{
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
	double sigma,   //channel noise standard deviation
	double x_d[N],  //x_dash  = noisy data symbol   
	double p1_d[N], //p1_dash = noisy parity#1 symbol
	double p2_d[N], //p2_dash = noisy parity#2 symbol
	double L_h[N],  //L_hat = likelihood of databit given entire observation
	int    X_h[N]   //X_hat = sliced MAP decisions
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
		sova(Lc, Le2_ip, x_d, p1_d, Le1,1);
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
		sova( Lc, Le1_p,  x_d_p, p2_d, Le2,0);
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


/*
gcc turbo_example.c -lm -o t; t
*/

void test(int *input, int *out, int *length)
{
	int    P1[N];     //encoder #1 parity bits
	int    P2[N];     //encoder #2 parity bits
	double x[N];      //databit mapped to symbol
	double p1[N];     //encoder #1 parity bit mapped to symbol
	double p2[N];     //encoder #2 parity bit mapped to symbol
	double x_d[N];    //x_dash  = noisy data symbol
	double p1_d[N];   //p1_dash = noisy parity#1 symbol
	double p2_d[N];   //p2_dash = noisy parity#2 symbol
	double L_h[N];    //L_hat = likelihood of databit given entire observation
	int    X_h[N];    //X_hat = sliced MAP decisions
	int    k;         //databit index (trellis stage)
	#if NOISEOFF == 1
	double noise[3*N];//noise from file
	FILE * pFile;	
	char c[51];
	#endif
	for(int i = 0;i < 8; i++){
		X[i] = input[i];
		printf("input[%i]= %i\n",i,input[i]);
	}
		

    /********************************
     *         INITIALISE           *
     ********************************/

	srand(1);    //init random number generator

    gen_tab();   //generate trellis tables


    /********************************
     *           ENCODER            *
     ********************************/

    turbo_encode(X, P1, P2);

      //map bits to symbols
	for(k = 0; k < N; k++) //for entire block length
	{
		x[k]  = X[k]  ? +1 : -1;  //map databit to symbol
		p1[k] = P1[k] ? +1 : -1;  //map parity #1 to symbol for(k = 0; k < N; k++)
		p2[k] = P2[k] ? +1 : -1;  //map parity #2 to symbol
	}
	#if DEBUG ==2
	for(k = 0; k < N; k++)
		printf("\nx=%f\tp1=%f\tp2=%f",x[k],p1[k],p2[k]);
	#endif

    /********************************
     *           CHANNEL            *
     ********************************/
	#if NOISEOFF == 1
	pFile = fopen (noise_file,"r");
	if (pFile!=NULL)
	{
		for(k=0; (k< 3*N) ;k++) 
			if (fgets(c, 50, pFile)!=NULL)
				noise[k] = atof(c);
			else
			{
				printf("problen with noise file %s", noise_file);
				//scanf("%c",c); 
				exit(1);
			}
	}
	else
	{
		printf("cant open file %s", noise_file);
		//scanf("%c",c);
		exit(1);
	}
	#endif
    //add some AWGN
	for(k = 0; k < N; k++)
	{
		#if NOISEOFF == 0
		x_d[k]  = x[k];
		p1_d[k] = p1[k];
		p2_d[k] = p2[k];
		#elif NOISEOFF == 1		
		x_d[k]  = x[k] + noise[k];
		p1_d[k] = p1[k] + noise[k+N];
		p2_d[k] = p2[k] + noise[k+2*N];
		#else
		x_d[k]  = x[k]  + sigma*normal();
		p1_d[k] = p1[k] + sigma*normal();
		p2_d[k] = p2[k] + sigma*normal();
		#endif
	}

    #if DEBUG ==2
	for(k = 0; k < N; k++)
		printf("\nx'=%f\tp1'=%f\tp2'=%f",x_d[k],p1_d[k],p2_d[k]);
	#endif
		

    /********************************
     *           DECODER            *
     ********************************/

    turbo_decode(1, x_d, p1_d, p2_d, L_h, X_h);

    //print soft decisions
	for(k = 0; k < N ; k++)
		printf("L_h[%i] = %f\n", k, L_h[k]);
	printf("\n");

    //print hard decisions
	printf("X_h = ");
    for(k = 0; k < N; k++)
    	printf("%i", X_h[k]);
    for(k = 0; k < N; k++)
    	*(out+k) = X_h[k];
	printf("\n");
}
// end of file: turbo_example.c

