#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>

#define DEBUG 1

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

void printMatrix(int rows, int cols, int matrix[][cols]) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			printf("%i ",matrix[i][j]);
		}
		printf("\n");
	}
}


void generateMatrices
(	// inputs
	int N,			// number of output symbols per input symbol (coderate r = 1/n)
	int M,			// constraint length, encoder memory length
	int* generator,	// generator polynoms, one for each output symbol
	// outputs
	int* vnextState,		// matrix of nextState in vector form (dimensions: NUM_STATES x 2)
	int* vpreviousState,	// 3D matrix of previousState in vector form (dimensions: 2 x NUM_STATES x 2)
	int* voutput			// matrix of output symbols in vector form (dimensions: NUM_STATES x 2)
) {
	
	const int NUM_STATES = pow(2,M);			// number of states: 2^M
	const int MAX_GENERATOR = pow(2,M+1) - 1;	// maximum number for a generator polynom
	// TODO: Check generators values for MAX_GENERATOR
	
	int nextState[NUM_STATES][2];			// nextState[current_state][input]
	int previousState[2][NUM_STATES][2];	// previousState[decision][current_state][past_input]
	int output[NUM_STATES][2];				// output[current_state][input]
	
	// initialization of previousState matrix
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < NUM_STATES; j++) {
			for (int k = 0; k < 2; k++) {
				previousState[i][j][k] = -1;
			}
		}
	}
	
	for (int state = 0; state < NUM_STATES; state++) {
		for (int input = 0; input < 2; input++) {
			// current_state represents the encoder state plus the input symbol as MSB! (encoder diagram view)
			int current_state = state | (input << M);
			
			// LSB of generator polynoms handles input symbol, therefore current_state is turned round!
			// maximum number of bits for current_state and the generator polynoms are M+1 (M state bits + 1 input bit)
			int turned_state = turnBitsRound(current_state, M+1);
			
			// calculate output for given state and input
			int out = 0;
		  	for (int i = 0; i < N; i++) {
		  		// temp is the output of symbol number i (defined by generator[i]) given current state turned_state
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
		for (int input = 0; input < 2; input++) {
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
	
	// bring matrices from matrix form to vector form
	for (int i = 0; i < NUM_STATES; i++) {
		for (int j = 0; j < 2; j++) {
			vnextState[2*i + j] = nextState[i][j]; 
			voutput[2*i + j] = output[i][j];  
		}
	}
	
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < NUM_STATES; j++) {
			for (int k = 0; k < 2; k++) {
				vpreviousState[NUM_STATES*2*i + 2*j + k] = previousState[i][j][k];
			}
		}	
	}
	
	#if DEBUG == 1
	printf("\nNext State\n");
	printMatrix(NUM_STATES, 2, nextState);
	printf("\nOutput\n");
	printMatrix(NUM_STATES, 2, output);
	printf("\nPrevious[0]\n");
	printMatrix(NUM_STATES, 2, previousState[0]);
	printf("\nPrevious[1]\n");
	printMatrix(NUM_STATES, 2, previousState[1]);
	#endif
}

/**********************
**      ENCODER      **
**********************/
void convolution_encode
(	// inputs
	int* input,
	int inputLen,
	int N,
	int M,
	// int nextState[][2],	// correct passing of a matrix
	int* vnextState,
	// int output[][2],
	int* voutput,
	// outputs
	int* code
) {
	const int NUM_STATES = pow(2,M);
	
	int nextState[NUM_STATES][2];			// nextState[current_state][input]
	int output[NUM_STATES][2];				// output[current_state][input]
	
	// bring matrices from vector form to matrix form
	for (int i = 0; i < NUM_STATES; i++) {
		for (int j = 0; j < 2; j++) {
			nextState[i][j] = vnextState[2*i + j]; 
			output[i][j] = voutput[2*i + j];  
		}
	}
	
	// convolution encoder
	int state = 0;
	int index = 0;
	for (int i = 0; i < inputLen; i++) {
		
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
	
	// bits are mapped: {0,1} --> {-1,+1}
	const int codeLen = N*(inputLen+M);
	for (int i = 0; i < codeLen; i++) {
		code[i] = 2 * code[i] - 1;
	}
	
	#if DEBUG == 1
	// print coded message
	for (int i = 0; i < codeLen; i++) {
		printf("%i",code[i]);
	}
	printf("\n");
	#endif
}

/**********************
**      DECODER      **
**********************/
void convolution_decode
(	// inputs
	int* code,
	int codeLen,
	int N,
	int M,
	int NUM_STATES,
	//int nextState[][2],
	int* vnextState,
	//int previousState[][NUM_STATES][2],	// correct passing of a matrix
	int* vpreviousState,
	//int output[][2],
	int* voutput,
	// outputs
	double* softOutput,
	int* hardOutput
) {

	const int msgLen = (codeLen / N);		// includes termination bits (M termination bits)
	
	int nextState[NUM_STATES][2];			// nextState[current_state][input]
	int previousState[2][NUM_STATES][2];	// previousState[decision][current_state][past_input]
	int output[NUM_STATES][2];				// output[current_state][input]
	
	double deltamin;
	double delta[msgLen+1][NUM_STATES];
	double metric[msgLen+1][NUM_STATES];
	int survivorBit[msgLen+1][NUM_STATES];
	int survivorStates[msgLen+1];
	int previous_matrix_decision[msgLen+1][NUM_STATES];
	
	// bring matrices from vector form to matrix form
	for (int i = 0; i < NUM_STATES; i++) {
		for (int j = 0; j < 2; j++) {
			nextState[i][j] = vnextState[2*i + j]; 
			output[i][j] = voutput[2*i + j];  
		}
	}
	// i = Matrixauswahl, j = State, k = Input-Bit
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < NUM_STATES; j++) {
			for (int k = 0; k < 2; k++) {
				previousState[i][j][k] = vpreviousState[NUM_STATES*2*i + 2*j + k];
			}
		}	
	}
	// NOTE: arrays have dimension msgLen+M+1 because i.e. delta[1][state] is delta from symbol at t=0 to t=1
	//		 so delta from symbol at t=msgLen+M-1 to t=msgLen+M (+M because of termination) is stored in
	//		 delta[msgLen+M][state] -> therefore delta[msgLen+M+1][NUM_STATES]
	//		 index 0 is not used!
	
	// initialization
	metric[0][0] = 0;
	for (int i = 1; i < NUM_STATES; i++) {
		metric[0][i] = -1000;
	}
	for (int i = 0; i < msgLen+1; i++) {
		for (int j = 0; j < NUM_STATES; j++) {
			previous_matrix_decision[i][j] = -1;
			survivorBit[i][j] = -1;	
		}
	}
	
	int Lc = 1;	// ACHTUNG!!!!!!
	int index = 0;	// index to select correct code bit
	
	// loop: time t
	for (int t = 1; t < msgLen+1; t++) {
		// loop: state s
		for (int s = 0; s < NUM_STATES; s++) {
			double Max[2];
			
			if (previousState[0][s][0] == -1 || previousState[0][s][1] == -1) {
				// Max[0] erhält die Metrik für previous State aus Matrix 0
				// Max[1] erhält die Metrik für previous State aus Matrix 1
				if (previousState[0][s][0] == -1) {
					// no previous for input bit 0
					survivorBit[t][s] = 1;
				}
				else if (previousState[0][s][1] == -1) {
					// no previous for input bit 1
					survivorBit[t][s] = 0;
				}
				Max[0] = metric[t-1][previousState[0][s][survivorBit[t][s]]];
				Max[1] = metric[t-1][previousState[1][s][survivorBit[t][s]]];
				for (int n = 0; n < N; n++) {
					int sr = N-n-1;	// shift right
					Max[0] += (((output[previousState[0][s][survivorBit[t][s]]][survivorBit[t][s]] >> sr) & 0x01) * 2 - 1) * Lc * code[index+n];
					Max[1] += (((output[previousState[1][s][survivorBit[t][s]]][survivorBit[t][s]] >> sr) & 0x01) * 2 - 1) * Lc * code[index+n];
				}
				
				previous_matrix_decision[t][s] = (Max[0] > Max[1]) ? 0 : 1;
			}
			else {
				// there are previous states for input bits 0 and 1
				// Max[0] erhält die Metrik für previous State mittels previous bit 0
				// Max[1] erhält die Metrik für previous State mittels previous bit 1
				for (int i = 0; i < 2; i++) {
					Max[i] = metric[t-1][previousState[0][s][i]];
					for (int n = 0; n < N; n++) {
						int sr = N-n-1;	// shift right
						Max[i] += (((output[previousState[0][s][i]][i] >> sr) & 0x01) * 2 - 1) * Lc * code[index+n];
					}
				}
				
				survivorBit[t][s] = (Max[0] > Max[1]) ? 0 : 1;
			}
			
			
			if (Max[0] > Max[1]) {
				metric[t][s] = Max[0];
				delta[t][s] = (Max[0] - Max[1]) / 2;
			}
			else {
				metric[t][s] = Max[1];
				delta[t][s] = (Max[1] - Max[0]) / 2;
			}
			
			#if DEBUG == 1
			printf("\ndelta[%d][%d]=%f\t M0=%.2f \tM1=%.2f \tmax=%f",t,s,delta[t][s], Max[0],Max[1], metric[t][s]);
			#endif
		}
		
		index += N;
	}
	
	int IsTerminated = 1;
	
	if (IsTerminated > 0) {
		survivorStates[msgLen] = 0;	
	}
	else {
		double max_metric = metric[msgLen][0];
		survivorStates[msgLen] = 0;
		for (int s = 1; s < NUM_STATES; s++) {
			if (max_metric < metric[msgLen][s]) {
				max_metric = metric[msgLen][s];
				survivorStates[msgLen] = s;
			}
		}
	}
	
	// reconstruction of survivor states, from right to left
	// survivor state at index msgLen is known, start loop at 2nd but last index (= msgLen-1)
	for (int t = msgLen-1; t >= 0; t--) {
		if (previous_matrix_decision[t+1][survivorStates[t+1]] == -1) {
			survivorStates[t] = previousState[0][survivorStates[t+1]][survivorBit[t+1][survivorStates[t+1]]];
		}
		else {
			survivorStates[t] = previousState[previous_matrix_decision[t+1][survivorStates[t+1]]][survivorStates[t+1]][survivorBit[t+1][survivorStates[t+1]]];
		}
	}
	
	int s;
	for (int t = M+1; t < msgLen+1; t++) {
		deltamin = delta[t][survivorStates[t]];
		// s ... state aus dem man in den survivor state kommen hätte können (2. Möglichkeit)
		
		if (previous_matrix_decision[t][survivorStates[t]] == -1) {
			//wenn keine Entscheidung getroffen wurde muss Bit flippen
			s = previousState[0][survivorStates[t]][survivorBit[t][survivorStates[t]] * (-1) + 1];
		}
		else {
			//falls Entscheidung getroffen wird muss State flippen
			s = previousState[previous_matrix_decision[t][survivorStates[t]]*(-1)+1][survivorStates[t]][survivorBit[t][survivorStates[t]]];
		}
		
		for (int i = t-1; i > 0; i--) {
			//sucht das kleinste delta auf dem survivor state
			if (delta[i][survivorStates[i]] < deltamin) {
				deltamin = delta[i][survivorStates[i]];
			}
			if (survivorBit[i][survivorStates[i]] != survivorBit[i][s]) {
				delta[i][survivorStates[i]] = deltamin;
				#if DEBUG == 1
				printf("\nupdate k=%d i=%d, delta=%f",t,i,delta[i][survivorStates[i]]);
				#endif
			}
			
			if (previous_matrix_decision[i][s] == -1) {
				//wenn keine Entscheidung getroffen wurde muss Bit flippen
				s = previousState[0][s][survivorBit[i][s]];
			} 
			else {
				//falls Entscheidung getroffen wird muss State flippen
				s = previousState[previous_matrix_decision[i][s]][s][survivorBit[i][s]];
			}
		}

	}

	for(int t = 1; t < msgLen+1; t++) {
		softOutput[t-1] = delta[t][survivorStates[t]] * (survivorBit[t][survivorStates[t]]*2-1);
		hardOutput[t-1] = (softOutput[t-1] > 0) ? 1 : 0;
	}
	
	#if DEBUG == 1
	printf("\n");
	for (int i = 0; i < msgLen; i++) {
		printf("%f\n",softOutput[i]);
	}
	for (int i = 0; i < msgLen; i++) {
		printf("%i",hardOutput[i]);
	}
	#endif
}

void main_temp() {
	int M = 2;
	int N = 2;
	int generator[2] = {7,5};
	
	int num_states = pow(2,M);
	
	int next[num_states*2];
	int prev[num_states*2*2];
	int out[num_states*2];
	
	generateMatrices(N,M,generator,next,prev,out);
	
	int msgLen = 12;
	int message[12] = {0,1,1,1,0,0,1,0,1,1,0,1};
	int codeLen = (msgLen + M) * N;
	int code[codeLen];
	
	convolution_encode(message, msgLen, N, M, next, out, code);
	
	code[4] = -code[4];
	code[14] = -code[14];
	code[15] = -code[15];
	code[20] = -code[20];
	
	double msgSoft[msgLen];
	int msgHard[msgLen];
	
	convolution_decode(code, codeLen, N, M, num_states, next, prev, out, msgSoft, msgHard);
	
	printf("\n");
	for (int i = 0; i < msgLen; i++) {
		printf("%i",message[i]);
	}
	printf("00");
}

void wrapper_generateMatrices
(	// inputs
	int* N,			// number of output symbols per input symbol (coderate r = 1/n)
	int* M,			// constraint length, encoder memory length
	int* generator,	// generator polynoms, one for each output symbol
	// outputs
	int* vnextState,		// matrix of nextState in vector form (dimensions: NUM_STATES x 2)
	int* vpreviousState,	// 3D matrix of previousState in vector form (dimensions: 2 x NUM_STATES x 2)
	int* voutput			// matrix of output symbols in vector form (dimensions: NUM_STATES x 2)
) {
	generateMatrices(*N, *M, generator, vnextState, vpreviousState, voutput);
}

void wrapper_convolution_encode
(	// inputs
	int* input,
	int* inputLen,
	int* N,
	int* M,
	// int nextState[][2],	// correct passing of a matrix
	int* vnextState,
	// int output[][2],
	int* voutput,
	// outputs
	int* code
) {
	convolution_encode(input, *inputLen, *N, *M, vnextState, voutput, code);
}

void wrapper_convolution_decode
(	// inputs
	int* code,
	int* codeLen,
	int* N,
	int* M,
	int* NUM_STATES,
	//int nextState[][2],
	int* vnextState,
	//int previousState[][NUM_STATES][2],	// correct passing of a matrix
	int* vpreviousState,
	//int output[][2],
	int* voutput,
	// outputs
	double* softOutput,
	int* hardOutput
) {
	convolution_decode(code, *codeLen, *N, *M, *NUM_STATES, vnextState, vpreviousState, voutput, softOutput, hardOutput);
}
