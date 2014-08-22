

/******************************************************************************************************************************************************************
***************************************************************    Author: Sai Prasanna.AN    *********************************************************************
***************************************************************      Date: July 15 2014       *********************************************************************
*******************************************************************************************************************************************************************  
*******************************************************************************************************************************************************************
************************************************************ This code is the Perceptron predictor ****************************************************************
*******************************************************************************************************************************************************************/


#ifndef PREDICTOR_H_SEEN
#define PREDICTOR_H_SEEN

#include <cstddef>
#include <cstdlib>
#include <bitset>
#include <inttypes.h>
#include "op_state.h"		
#include "tread.h"		
#include "math.h"		


#define ASSERT(cond) if (!(cond)) {printf("assert line %d\n",__LINE__); exit(EXIT_FAILURE);}


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  MACRO DEFINITION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


#define PERCEPTRON_HISTORY	62
#define NUM_PERCEPTRONS		2731
#define PERCEPTRON_BITS		8
#define MAX_WEIGHT		((1<<(PERCEPTRON_BITS-1))-1)
#define MIN_WEIGHT		(-(MAX_WEIGHT+1))
#define THETA			((int) (1.93 * PERCEPTRON_HISTORY + 14))
#define NUM_UPDATE_ENTRIES	100

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  END OF MACROS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


using namespace std;

typedef uint32_t address_t;
typedef bitset <MAXHIST> history_t;


/* This class defines cyclic shift register for folding long global history bits to smaller number of bits   * 
 * Original length of history bits for TAGE tables form a geometric series                                   *
 * Compressed length if used for indexing will then be LOGG number of bits                                   *
 * Compressed length if used for tag match, depends on number of tag bits for a given table                  */



/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  ALL THE PREDICTOR IS HERE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

class PREDICTOR
{
  public:



/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  PREDICTOR STRORAGE DATA - PERCEPTRON  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  typedef struct {
	int weights[PERCEPTRON_HISTORY+1];
  } perceptron;

  typedef struct {
	char dummy_counter;
	int prediction, output;
	unsigned long long int history;
	perceptron *perc;
  } perceptron_state;

  perceptron perceptrons[NUM_PERCEPTRONS];
  perceptron_state perceptron_state_buf[NUM_UPDATE_ENTRIES];
  int perceptron_state_buf_ctr;
  unsigned long long int spec_global_history, global_history;
  perceptron_state temp;
  perceptron_state *u;
   
 
 /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  CONSTRUCTOR INITIALIZE DATA  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
 
   void initialize_perceptron (perceptron *p) {
         for (int i=0; i<=PERCEPTRON_HISTORY; i++) p->weights[i] = 0;
     }

     void initialize_perceptron_predictor (void) {
          spec_global_history = 0;
          global_history = 0;
          perceptron_state_buf_ctr = 0;
          for (int i=0; i<NUM_PERCEPTRONS; i++)
             initialize_perceptron (&perceptrons[i]);
     }

 PREDICTOR ()
 {
    
    ghist = 0;
    spec_global_history = 0;
    global_history = 0;
    perceptron_state_buf_ctr = 0;
    initialize_perceptron_predictor();


 }

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  END OF CONSTRUCTOR >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/


 
 
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  getperceptron  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  /* Prediction by the perceptron predictor */
  
  bool getperceptron(address_t pc)
  {
        uint32_t address = pc;     /* Address of the corresponding branch instruction */
    	int index, i, output, *w;
	unsigned long long int mask;
	perceptron *p;

	/* get a pointer to the next "free" perceptron_state,
	* bumping up the pointer (and possibly letting it wrap around) 
	*/

	u = &perceptron_state_buf[perceptron_state_buf_ctr++];
	if (perceptron_state_buf_ctr >= NUM_UPDATE_ENTRIES)
		perceptron_state_buf_ctr = 0;

	/* hash the address to get an index into the table of perceptrons */

	index = address % NUM_PERCEPTRONS;

	/* get pointers to that perceptron and its weights */

	p = &perceptrons[index];
	w = &p->weights[0];

	/* initialize the output to the bias weight, and bump the pointer
	* to the weights
	*/

	output = *w++;

	/* find the (rest of the) dot product of the history register and the perceptron weights.  note that, instead of actually
	* doing the expensive multiplies, we simply add a weight when the corresponding branch in the history register is taken, or
	* subtract a weight when the branch is not taken.  this also lets us use binary instead of bipolar logic to represent the history register
	*/

	for (mask=1,i=0; i<PERCEPTRON_HISTORY; i++,mask<<=1,w++) {
	if (spec_global_history & mask)
	output += *w;
	else
	output += -*w;
	}

	/* record the various values needed to update the predictor */

	u->output = output;
	u->perc = p;
	u->history = spec_global_history;
	u->prediction = output >= 0;
	u->dummy_counter = u->prediction ? 3 : 0;

	/* update the speculative global history register */
	
	spec_global_history <<= 1;
	spec_global_history |= u->prediction;

        return u->dummy_counter;



  } 

  
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  perceptronupdate  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* Base perceptron predictor update */   

  void perceptronupdate (address_t pc, bool taken)
  {
        int i, y, *w;
	unsigned long long int mask, history;

	/* update the real global history shift register */

	global_history <<= 1;
	global_history |= taken;

	/* if this branch was mispredicted, restore the speculative history to the last known real history */

	if (u->prediction != taken) spec_global_history = global_history;


	/* if the output of the perceptron predictor is outside of the range [-THETA,THETA] *and* the prediction was correct, 
	 * then we don't need to adjust the weights */

	if (u->output > THETA)
	  y = 1;
	else if (u->output < -THETA)
	  y = 0;
	else
	  y = 2;
	if (y == 1 && taken) return;
	if (y == 0 && !taken) return;

	/* w is a pointer to the first weight (the bias weight) */

	w = &u->perc->weights[0];

	/* if the branch was taken, increment the bias weight, else decrement it, with saturating arithmetic */

	if (taken)
	  (*w)++;
	else
	  (*w)--;
	
	if (*w > MAX_WEIGHT) *w = MAX_WEIGHT;
	if (*w < MIN_WEIGHT) *w = MIN_WEIGHT;

 	/* now w points to the next weight */

	w++;

	/* get the history that led to this prediction */

	history = u->history;

	/* for each weight and corresponding bit in the history register... */

	for (mask=1,i=0; i<PERCEPTRON_HISTORY; i++,mask<<=1,w++) {

	  /* if the i'th bit in the history positively correlates with this branch outcome, increment the corresponding 
           * weight, else decrement it, with saturating arithmetic */

	  if (!!(history & mask) == taken) {
	    (*w)++;
	  if (*w > MAX_WEIGHT) *w = MAX_WEIGHT;
	  } else {
	    (*w)--;
	  if (*w < MIN_WEIGHT) *w = MIN_WEIGHT;
	  }
	
        }

  }  

};
#endif // PREDICTOR_H_SEEN

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  END ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
