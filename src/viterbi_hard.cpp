#include <Rcpp.h>
using namespace Rcpp;

#define DEBUG 1

// [[Rcpp::export]]
IntegerVector c_convolutionDecode_hard(IntegerVector code, int N, int M, IntegerMatrix previousState, IntegerMatrix output) {
	
	const int codeLen = code.size();
	const int msgLen = (codeLen / N);		// includes termination bits (M termination bits)
	const int NUM_STATES = pow(2,M);
	
	int metric[msgLen+1][NUM_STATES];
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
		metric[0][i] = 10;
	}
	for (int i = 0; i < msgLen+1; i++) {
		for (int j = 0; j < NUM_STATES; j++) {
			previousMatrixColumn[i][j] = -1;
		}
	}
	
	// reverse mapping for hard decision decoding {-1,+1} --> {0,1}
	for (int i = 0; i < codeLen; i++) {
		code[i] = (code[i] + 1) / 2;
	}
	
	int index = 0;	// index to select correct code bit
	// metric calculation
	// loop: time t
	for (int t = 1; t < msgLen+1; t++) {
		// loop: state s
		for (int s = 0; s < NUM_STATES; s++) {
			int Min[2];
			
			if (previousState(s,0) == -1 || previousState(s,1) == -1) {
				if (previousState(s,0) == -1) {
					// no previous for input bit 0
					survivorBit[t][s] = 1;
				}
				else if (previousState(s,1) == -1) {
					// no previous for input bit 1
					survivorBit[t][s] = 0;
				}
				
				Min[0] = metric[t-1][previousState(s,survivorBit[t][s])];
				Min[1] = metric[t-1][previousState(s,2)];
				for (int n = 0; n < N; n++) {
					int sr = N-n-1;	// shift right
					Min[0] += ((output(previousState(s,survivorBit[t][s]),survivorBit[t][s]) >> sr) & 0x01) ^ code[index+n];
					Min[1] += ((output(previousState(s,2),survivorBit[t][s]) >> sr) & 0x01) ^ code[index+n];
				}
				
				previousMatrixColumn[t][s] = (Min[0] < Min[1]) ? survivorBit[t][s] : 2;
			}
			else {
				// there are previous states for input bits 0 and 1
				for (int i = 0; i < 2; i++) {
					Min[i] = metric[t-1][previousState(s,i)];
					for (int n = 0; n < N; n++) {
						int sr = N-n-1;	// shift right
						Min[i] += ((output(previousState(s,i),i) >> sr) & 0x01) ^ code[index+n];
					}
				}
				
				survivorBit[t][s] = (Min[0] < Min[1]) ? 0 : 1;
			}
			
			metric[t][s] = (Min[0] < Min[1]) ? Min[0] : Min[1];
			
			#if DEBUG == 1
			printf("metric[%d][%d]=%d\t M0=%d \tM1=%d\n",t,s,metric[t][s], Min[0],Min[1]);
			#endif
		}
		
		index += N;
	}
	
	#if DEBUG == 1
	for (int j = 0; j < NUM_STATES; j++) {
		for (int i = 0; i < msgLen + 1; i++) {
			printf("%3d",metric[i][j]);
		}
		printf("\n");
	}
	#endif
	
	int IsTerminated = 1;
	
	if (IsTerminated > 0) {
		survivorStates[msgLen] = 0;	
	}
	else {
		int min_metric = metric[msgLen][0];
		survivorStates[msgLen] = 0;
		for (int s = 1; s < NUM_STATES; s++) {
			if (min_metric > metric[msgLen][s]) {
				min_metric = metric[msgLen][s];
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
			survivorStates[t] = previousState(survivorStates[t+1], previousMatrixColumn[t+1][survivorStates[t+1]]);
		}
	}
	
	IntegerVector hardOutput(msgLen);
	
	for(int t = 1; t < msgLen+1; t++) {
		hardOutput[t-1] = survivorBit[t][survivorStates[t]];
	}
	
	#if DEBUG == 1
	for (int i = 0; i < msgLen; i++) {
		printf("%i",hardOutput[i]);
	}
	#endif
	
	return hardOutput;
}
