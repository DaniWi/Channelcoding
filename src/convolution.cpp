#include <Rcpp.h>
using namespace Rcpp;

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

// [[Rcpp::export]]
List c_generateMatrices(int N, int M, IntegerVector generator) {

	const int NUM_STATES = pow(2,M);			// number of states: 2^M
	
	IntegerMatrix nextState(NUM_STATES,2);
	IntegerMatrix output(NUM_STATES,2);
	IntegerMatrix previousState(NUM_STATES,3);
	
	// initialization of previousState matrix
	for (int i = 0; i < NUM_STATES; i++) {
		for (int j = 0; j < 2; j++) {
			previousState(i,j) = -1;
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
	    		// cast to int of generator, because generator[i] is a Numeric
	    		int temp = sumDigits(turned_state & generator[i], 2) % 2;
	    		out = (out << 1) | temp;
	  		}
			output(state, input) = out;
		
			// calculate new state
			nextState(state, input) = current_state >> 1;
		}
	}
	
	// filling of previousState matrix
	for (int state = 0; state < NUM_STATES; state++) {
		for (int input = 0; input < 2; input++) {
			if (previousState(nextState(state, input), input) == -1) {
				// 1st matrix's spot is still free
				previousState(nextState(state, input), input) = state;
			}
			else {
				// there are two states that have the same input bit for transition to state
				// --> use 2nd spot to save other previous state with same input bit
				previousState(nextState(state, input), 2) = state;
			}
		}
	}
	
	List matrices = List::create(Rcpp::Named("nextState") = nextState,
								 Rcpp::Named("prevState") = previousState,
								 Rcpp::Named("output") = output);

	return matrices;
}

// [[Rcpp::export]]
IntegerVector c_convolutionEncode(IntegerVector input, int N, int M, IntegerMatrix nextState, IntegerMatrix output) {
	
	const int inputLen = input.size();
	const int codeLen = (inputLen + M) * N;
	
	IntegerVector code(codeLen);
	
	int state = 0;
	int index = 0;
	for (int i = 0; i < inputLen; i++) {
		
		int out = output(state, input[i]);
		
		#if DEBUG == 1
		printf("state %i, input %i, output %i%i, new state %i\n",state,input[i],(out>>1)&1,out&1,nextState(state, input[i]));
		#endif
		
		state = nextState(state, input[i]);
		// output consists of N bits per coded symbol
		for (int j = N-1; j >= 0; j--) {
			code[index] = (out >> j) & 0x01;
			index++;
		}
	}
	
	// termination: input of M times 0
	for (int i = 0; i < M; i++) {
		int out = output(state,0);
		state = nextState(state,0);
		for (int j = N-1; j >= 0; j--) {
			code[index] = (out >> j) & 0x01;
			index++;
		}
	}
	
	// bits are mapped: {0,1} --> {-1,+1}
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
	
	return code;
}

// [[Rcpp::export]]
List c_convolutionDecode(NumericVector code, int N, int M, IntegerMatrix previousState, IntegerMatrix output) {
	
	const int codeLen = code.size();
	const int msgLen = (codeLen / N);		// includes termination bits (M termination bits)
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
	for (int i = 1; i < NUM_STATES; i++) {
		metric[0][i] = -1000;
	}
	for (int i = 0; i < msgLen+1; i++) {
		for (int j = 0; j < NUM_STATES; j++) {
			previousMatrixColumn[i][j] = -1;
		}
	}
	
	int Lc = 1;	// ACHTUNG!!!!!!
	int index = 0;	// index to select correct code bit
	
	// loop: time t
	for (int t = 1; t < msgLen+1; t++) {
		// loop: state s
		for (int s = 0; s < NUM_STATES; s++) {
			double Max[2];
			
			if (previousState(s,0) == -1 || previousState(s,1) == -1) {
				// Max[0] erh�lt die Metrik f�r previous State aus Spalte 0 bzw 1
				// Max[1] erh�lt die Metrik f�r previous State aus Spalte 2
				if (previousState(s,0) == -1) {
					// no previous for input bit 0
					survivorBit[t][s] = 1;
				}
				else if (previousState(s,1) == -1) {
					// no previous for input bit 1
					survivorBit[t][s] = 0;
				}
				Max[0] = metric[t-1][previousState(s,survivorBit[t][s])];
				Max[1] = metric[t-1][previousState(s,2)];
				for (int n = 0; n < N; n++) {
					int sr = N-n-1;	// shift right
					Max[0] += (((output(previousState(s,survivorBit[t][s]),survivorBit[t][s]) >> sr) & 0x01) * 2 - 1) * Lc * code[index+n];
					Max[1] += (((output(previousState(s,2),survivorBit[t][s]) >> sr) & 0x01) * 2 - 1) * Lc * code[index+n];
				}
				
				previousMatrixColumn[t][s] = (Max[0] > Max[1]) ? survivorBit[t][s] : 2;
			}
			else {
				// there are previous states for input bits 0 and 1
				// Max[0] erh�lt die Metrik f�r previous State mittels previous bit 0
				// Max[1] erh�lt die Metrik f�r previous State mittels previous bit 1
				for (int i = 0; i < 2; i++) {
					Max[i] = metric[t-1][previousState(s,i)];
					for (int n = 0; n < N; n++) {
						int sr = N-n-1;	// shift right
						Max[i] += (((output(previousState(s,i),i) >> sr) & 0x01) * 2 - 1) * Lc * code[index+n];
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
		if (previousMatrixColumn[t+1][survivorStates[t+1]] == -1) {
			survivorStates[t] = previousState(survivorStates[t+1], survivorBit[t+1][survivorStates[t+1]]);
		}
		else {
			//survivorStates[t] = previousState[previousMatrixColumn[t+1][survivorStates[t+1]]][survivorStates[t+1]][survivorBit[t+1][survivorStates[t+1]]];
			survivorStates[t] = previousState(survivorStates[t+1], previousMatrixColumn[t+1][survivorStates[t+1]]);
		}
	}
	
	int s;
	for (int t = M+1; t < msgLen+1; t++) {
		deltamin = delta[t][survivorStates[t]];
		// s ... state aus dem man in den survivor state kommen h�tte k�nnen (2. M�glichkeit)
		
		if (previousMatrixColumn[t][survivorStates[t]] == -1) {
			// wenn keine Entscheidung getroffen wurde muss Bit flippen
			s = previousState(survivorStates[t], survivorBit[t][survivorStates[t]] * (-1) + 1);
		}
		else {
			// falls Entscheidung getroffen wird muss Spalte wechseln
			// 2 F�lle:	a) 0 <--> 2
			//			b) 1 <--> 2
			// ob 0 oder 1 steht im survivorBit
			int column;
			if (previousMatrixColumn[t][survivorStates[t]] == 2) {
				column = survivorBit[t][survivorStates[t]];
			}
			else {
				column = 2;
			}
			//s = previousState[previousMatrixColumn[t][survivorStates[t]]*(-1)+1][survivorStates[t]][survivorBit[t][survivorStates[t]]];
			s = previousState(survivorStates[t], column);
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
			
			if (previousMatrixColumn[i][s] == -1) {
				//wenn keine Entscheidung getroffen wurde muss Bit flippen
				s = previousState(s, survivorBit[i][s]);
			} 
			else {
				//falls Entscheidung getroffen wird muss Spalte wechseln
				s = previousState(s, previousMatrixColumn[i][s]);
			}
		}

	}

	NumericVector softOutput(msgLen);
	IntegerVector hardOutput(msgLen);
	
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
	
	List result = List::create(Rcpp::Named("softOutput") = softOutput,
							   Rcpp::Named("hardOutput") = hardOutput);
	
	return result;
}
