#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>

#define DEBUG 0

int sumDigits(unsigned int n, int base) {
    int sum = 0;
    for (; n; n /= base){
    	sum += n % base;
    }
    return sum;
}

int turnBitsRound(int num, int nbits) {
	int result = 0;
	for (int i=0; i < nbits; i++) {
		result |= ((num >> i) & 0x01) << (nbits - i - 1);
	}
	return result;
}

void printMatrix(int matrix[][2], int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 2; j++) {
			printf("%i ",matrix[i][j]);
		}
		printf("\n");
	}
}

// outputs:
//		array: generator polynoms
void generateMatrices
(	int K,		// number of input symbols. K = 2 ... {0,1}
	int N,		// number of output symbols per input symbol (coderate r = 1/n)
	int M		// constraint length, encoder memory length
) {
	
	const int NUM_STATES = pow(2,M);			// number of states: 2^M
	const int MAX_GENERATOR = pow(2,M+1) - 1;	// maximum number for a generator polynom
	
	int generator[N];	// generator polynom for each output symbol
	for (int i = 0; i < N; i++) {
		printf("generator %i: ", i+1);
		scanf("%i",generator+i);
		generator[i] &= MAX_GENERATOR;	// if input is too high just consider bits within maximum generator
	}
	
	int previousState[K][NUM_STATES][K];	// previousState[decision][current_state][past_input]
	int nextState[NUM_STATES][K];			// nextState[current_state][input]
	int output[NUM_STATES][K];				// output[current_state][input]
	
	// initialization of previousState matrix
	for (int i = 0; i < K; i++) {
		for (int j = 0; j < NUM_STATES; j++) {
			for (int k = 0; k < K; k++) {
				previousState[i][j][k] = -1;
			}
		}
	}
	
	for (int state = 0; state < NUM_STATES; state++) {
		for (int input = 0; input < K; input++) {
			// current_state represents the encoder state plus the input symbol as MSB! (encoder diagram view)
			int current_state = state | (input << M);
			
			// LSB of generator polynoms handles input symbol, therefore current_state is turned round!
			// maximum number of bits for current_state and the generator polynoms are M+1 (M state bits + 1 input bit)
			int turned_state = turnBitsRound(current_state, M+1);
			
			// calculate output for given state and input
			int out = 0;
		  	for (int i = 0; i < N; i++) {
		  		int temp = sumDigits(turned_state & generator[i], 2) % 2;
		  		out = (out << 1) | temp;
			}
			output[state][input] = out;
			
			// calculate new state
			nextState[state][input] = current_state >> 1;
		}
	}
	
	// filling of previousState matrix
	for (int state = 0; state < NUM_STATES; state++) {
		for (int input = 0; input < K; input++) {
			if (previousState[0][nextState[state][input]][input] == -1) {
				// 1st matrix's spot is still free
				previousState[0][nextState[state][input]][input] = state;
			}
			else {
				// there are two states that have the same input bit for transition to state
				// --> use 2nd spot to save other previous state with same input bit
				previousState[1][nextState[state][input]][input] = state;
			}
			
		}
	}
	
	#if DEBUG == 1
	printf("\nNext State\n");
	printMatrix(nextState,NUM_STATES);
	printf("\nOutput\n");
	printMatrix(output,NUM_STATES);
	printf("\nPrevious\n");
	printMatrix(previousState,NUM_STATES);
	#endif
}


void convolution_encode
(	int* input,
	int input_len,
	int K,
	int N,
	int M,
	int** nextState,
) {
	
	const int code_len = N*(input_len+M);
	int code[code_len];
	
	// convolution encoder
	int state = 0;
	int index = 0;
	for (int i = 0; i < input_len; i++) {
		
		int out = output[state][input[i]];
		
		#if DEBUG == 1
		printf("state %i, input %i, output %i%i, new state %i\n",state,input[i],(out>>1)&1,out&1,nextState[state][input[i]]);
		#endif
		
		state = nextState[state][input[i]];
		// output consists of N bits per coded symbol
		for (int j = N-1; j >= 0; j--) {
			code[index] = (out >> j) & 0x01;
			index++;
		}
	}
	
	// termination: input of M times 0
	for (int i = 0; i < M; i++) {
		int out = output[state][0];
		state = nextState[state][0];
		for (int j = N-1; j >= 0; j--) {
			code[index] = (out >> j) & 0x01;
			index++;
		}
	}
	
	#if DEBUG == 1
	// print coded message
	for (int i = 0; i < code_len; i++) {
		printf("%i",code[i]);
	}
	printf("\n");
	#endif
}

void convolution_decode
(	int* input,
	int input_len,
	int K,
	int N,
	int M,
	int** nextState,
) {

	double metric[msg_len+M+1][NUM_STATES];
	double delta[msg_len+M+1][NUM_STATES];
	double deltamin;
	int survivor_bit[msg_len+M+1][NUM_STATES];
	int survivor_states[msg_len+M+1];
	
	// NOTE: arrays have dimension msg_len+M+1 because i.e. delta[1][state] is delta from symbol at t=0 to t=1
	//		 so delta from symbol at t=msg_len+M-1 to t=msg_len+M (+M because of termination) is stored in
	//		 delta[msg_len+M][state] -> therefore delta[msg_len+M+1][NUM_STATES]
	//		 index 0 is not used!
	
	// initialization
	metric[0][0] = 0;
	for (int i = 1; i < NUM_STATES; i++) {
		metric[0][i] = -1000;
	}
	
	index = 0;
	// int trellis_len = msg_len + M + 1;
	
	// ACHTUNG
	int Lc = 1;	// ACHTUNG!!!!!!
	// ACHTUNG
	
	// loop: time t
	for (int t = 0; t < msg_len+M+1; t++) {
		// loop: state s
		for (int s = 0; s < NUM_STATES; s++) {
			
			double Max[K];
			// loop: input i {0,1}
			for (int i = 0; i < K; i++) {
				
				Max[i] = metric[t-1][previousState[s][t]];
				// loop: 1 input symbol is coded to N output symbols
				// scalar product: N factors! sum the together!
				for (int n = 0; n < N; n++) {
					int sr = N-n-1;	// shift right
					Max[i] += (((output[previousState[s][t]][t] >> sr) & 0x01) * 2 - 1) * Lc * code[index+n-1];
					// La missing in previous calculation of Max[i]!!!
				}
				
				if (Max[0] > Max[1]) {
					survivor_bit[t][s] = 0;
					metric[t][s] = Max[0];
					delta[t][s] = (Max[0] - Max[1]) / 2;
				}
				else {
					survivor_bit[t][s] = 1;
					metric[t][s] = Max[1];
					delta[t][s] = (Max[1] - Max[0]) / 2;
				}
			}
		}
	}
	
	int IsTerminated = 1;
	
	if (IsTerminated > 0) {
		survivor_states[msg_len+M] = 0;	
	}
	else {
		double max_metric = metric[msg_len+M][0];
		survivor_states[msg_len+M] = 0;
		for(int s = 1; s < NUM_STATES; s++) {
			if (max_metric < metric[msg_len+M][s]) {
				max_metric = metric[msg_len+M][s];
				survivor_states[msg_len+M] = s;
			}
		}
	}
	
	// survivor state at index msg_len+M is known, start loop at 2nd but last index (= msg_len+M-1)
	for (int t = msg_len+M-1; t >= 0; t--) {
		survivor_states[t] = previousState[survivor_states[t+1]][survivor_bit[t+1][survivor_states[t+1]]];
	}
	for (int t = M+1; t < msg_len+M+1; t++) {
		deltamin = delta[t][survivor_states[t]];
		// s ... state aus dem man in den survivor state kommen hätte können (2. Möglichkeit)
		int s = previousState[ survivor_states[t] ][survivor_bit[t][survivor_states[t]] * (-1) + 1];
		for (int i = t-1; i > 0; i--) {
			//sucht das kleinste delta auf dem survivor state
			if (delta[i][survivor_states[i]] < deltamin) {
				deltamin = delta[i][survivor_states[i]];
			}
			if (survivor_bit[i][survivor_states[i]] != survivor_bit[i][t]) {
					delta[i][survivor_states[i]] = deltamin;
			}
			s = previousState[t][survivor_bit[i][t]];
		}

	}

	double Le[msg_len+M];
	for(int t = 1; t < msg_len+M+1; t++) {
		Le[t-1] = delta[t][survivor_states[t]] * (survivor_bit[t][survivor_states[t]]*2-1);
	}	
}

void main() {
	int K = 2;	// number of input symbols. K = 2 => {0,1}
	int N;		// number of output symbols per input symbol (coderate r = 1/n)
	int M;		// constraint length, encoder memory length
	
	// input of N
	do {
		printf("N: ");
		scanf("%i",&N);
	} while (N <= 0);
	
	// input of M
	do {
		printf("M: ");
		scanf("%i",&M);
	} while (M <= 0);
	
	const int NUM_STATES = pow(2,M);			// number of states: 2^M
	const int MAX_GENERATOR = pow(2,M+1) - 1;	// maximum number for a generator polynom
	printf("maximum generator polynom: %i\n", MAX_GENERATOR);
	
	int generator[N];	// generator polynom for each output symbol
	for (int i = 0; i < N; i++) {
		printf("generator %i: ", i+1);
		scanf("%i",generator+i);
		generator[i] &= MAX_GENERATOR;	// if input is too high just consider bits within maximum generator
	}
	
	int previousState[NUM_STATES][K];	// previousState[current_state][past_input]
	int nextState[NUM_STATES][K];		// nextState[current_state][input]
	int output[NUM_STATES][K];			// output[current_state][input]
	
	for (int state = 0; state < NUM_STATES; state++) {
		for (int input = 0; input < K; input++) {
			// current_state represents the encoder state plus the input symbol as MSB! (encoder diagram view)
			int current_state = state | (input << M);
			
			// LSB of generator polynoms handles input symbol, therefore current_state is turned round!
			// maximum number of bits for current_state and the generator polynoms are M+1 (M state bits + 1 input bit)
			int turned_state = turnBitsRound(current_state, M+1);
			
			// calculate output for given state and input
			int out = 0;
		  	for (int i = 0; i < N; i++) {
		  		int temp = sumDigits(turned_state & generator[i], 2) % 2;
		  		out = (out << 1) | temp;
			}
			output[state][input] = out;
			
			// calculate new state
			nextState[state][input] = current_state >> 1;
		}
	}
	
	for (int state = 0; state < NUM_STATES; state++) {
		for (int input = 0; input < K; input++) {
			previousState[nextState[state][input]][input] = state;
		}
	}
	
	
	printf("\nNext State\n");
	printMatrix(nextState,NUM_STATES);
	printf("\nOutput\n");
	printMatrix(output,NUM_STATES);
	printf("\nPrevious\n");
	printMatrix(previousState,NUM_STATES);
	
	
	int msg_len = 3;
	int code_len = N*(msg_len+M);
	int input[3] = {1, 0, 1};
	int code[code_len];
	
	// convolution encoder
	int state = 0;
	int index = 0;
	for (int i = 0; i < msg_len; i++) {
		
		int out = output[state][input[i]];
		printf("state %i, input %i, output %i%i, new state %i\n",state,input[i],(out>>1)&1,out&1,nextState[state][input[i]]);
		state = nextState[state][input[i]];
		// output consists of N bits per coded symbol
		for (int j = N-1; j >= 0; j--) {
			code[index] = (out >> j) & 0x01;
			index++;
		}
	}
	
	// termination: input of M times 0
	for (int i = 0; i < M; i++) {
		int out = output[state][0];
		state = nextState[state][0];
		for (int j = N-1; j >= 0; j--) {
			code[index] = (out >> j) & 0x01;
			index++;
		}
	}
	
	// print coded message
	for (int i = 0; i < code_len; i++) {
		printf("%i",code[i]);
	}
	printf("\n");
	
	
	// convolution decoder: viterbi algorithm
	double metric[msg_len+M+1][NUM_STATES];
	double delta[msg_len+M+1][NUM_STATES];
	double deltamin;
	int survivor_bit[msg_len+M+1][NUM_STATES];
	int survivor_states[msg_len+M+1];
	
	// NOTE: arrays have dimension msg_len+M+1 because i.e. delta[1][state] is delta from symbol at t=0 to t=1
	//		 so delta from symbol at t=msg_len+M-1 to t=msg_len+M (+M because of termination) is stored in
	//		 delta[msg_len+M][state] -> therefore delta[msg_len+M+1][NUM_STATES]
	//		 index 0 is not used!
	
	// initialization
	metric[0][0] = 0;
	for (int i = 1; i < NUM_STATES; i++) {
		metric[0][i] = -1000;
	}
	
	index = 0;
	// int trellis_len = msg_len + M + 1;
	
	// Channel reliability
	int Lc = 1;	// ACHTUNG!!!!!!
	
	// loop: time t
	for (int t = 0; t < msg_len+M+1; t++) {
		// loop: state s
		for (int s = 0; s < NUM_STATES; s++) {
			
			double Max[K];
			// loop: input i {0,1}
			for (int i = 0; i < K; i++) {
				
				Max[i] = metric[t-1][previousState[s][t]];
				// loop: 1 input symbol is coded to N output symbols
				// scalar product: N factors! sum the together!
				for (int n = 0; n < N; n++) {
					int sr = N-n-1;	// shift right
					Max[i] += (((output[previousState[s][t]][t] >> sr) & 0x01) * 2 - 1) * code[index+n-1] * Lc;
					// La missing in previous calculation of Max[i]!!!
				}
			}
			
			if (Max[0] > Max[1]) {
				survivor_bit[t][s] = 0;
				metric[t][s] = Max[0];
				delta[t][s] = (Max[0] - Max[1]) / 2;
			}
			else {
				survivor_bit[t][s] = 1;
				metric[t][s] = Max[1];
				delta[t][s] = (Max[1] - Max[0]) / 2;
			}
		}
	}
	
	int IsTerminated = 1;
	
	if (IsTerminated > 0) {
		survivor_states[msg_len+M] = 0;	
	}
	else {
		double max_metric = metric[msg_len+M][0];
		survivor_states[msg_len+M] = 0;
		for(int s = 1; s < NUM_STATES; s++) {
			if (max_metric < metric[msg_len+M][s]) {
				max_metric = metric[msg_len+M][s];
				survivor_states[msg_len+M] = s;
			}
		}
	}
	
	// survivor state at index msg_len+M is known, start loop at 2nd but last index (= msg_len+M-1)
	for (int t = msg_len+M-1; t >= 0; t--) {
		survivor_states[t] = previousState[survivor_states[t+1]][survivor_bit[t+1][survivor_states[t+1]]];
	}
	for (int t = M+1; t < msg_len+M+1; t++) {
		deltamin = delta[t][survivor_states[t]];
		// s ... state aus dem man in den survivor state kommen hätte können (2. Möglichkeit)
		int s = previousState[ survivor_states[t] ][survivor_bit[t][survivor_states[t]] * (-1) + 1];
		for (int i = t-1; i > 0; i--) {
			//sucht das kleinste delta auf dem survivor state
			if (delta[i][survivor_states[i]] < deltamin) {
				deltamin = delta[i][survivor_states[i]];
			}
			if (survivor_bit[i][survivor_states[i]] != survivor_bit[i][t]) {
					delta[i][survivor_states[i]] = deltamin;
			}
			s = previousState[t][survivor_bit[i][t]];
		}

	}

	double Le[msg_len+M];
	for(int t = 1; t < msg_len+M+1; t++) {
		Le[t-1] = delta[t][survivor_states[t]] * (survivor_bit[t][survivor_states[t]]*2-1);
	}
	
}
