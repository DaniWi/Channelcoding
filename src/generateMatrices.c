#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

void main() {
	const int K = 2;	// number of input symbols. K = 2 => {0,1}
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
	
}
