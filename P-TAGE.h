

/******************************************************************************************************************************************************************
***************************************************************    Author: Sai Prasanna.AN    *********************************************************************
***************************************************************      Date: July 20 2014       *********************************************************************
******************************************************************************************************************************************************************* 
*******************    The structure of the code is essentially derived  from the tagged PPM predictor simulator and the OGEHL predictor simulator   ************** 
*******************************************************************************************************************************************************************
************************************ This code is the P-TAGE predictor with 5/8/14 components based on the Macro definition ***************************************/


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


/* The predictor features NHIST tagged components + a base bimodal component                     *                                               
 * A tagged entry in tagged table Ti (0<= i< NHIST)  features a TBITS-(i+ (NHIST & 1))/2 tag,    *
 * a CBITS prediction counter and a 2-bit useful counter                                         *
 * Tagged components feature 2^LOGG entries                                                      *
 * On the bimodal table: hysteresis is shared among 4 counters                                   */


/* Prediction counter bits for TAGE tables */
#ifndef CBITS
#define CBITS 3
#endif


/* The default predictor:                                        *
 * Features 7 tagged components and a base bimodal component     *
 * NHIST = 7, LOGB =13, LOGG=9, CBITS=3                          */ 

#ifndef LOGB
#define LOGB 17
#endif

#ifndef NHIST
#define NHIST 7
#endif

#ifndef LOGG
#define LOGG (LOGB-4)
#endif


/* Different TAGE tables have different tag lengths
 * Total width of an entry in the tagged table with the longest history length   */ 

#ifndef TBITS
#define TBITS 12
#endif


/* We use Geometric history length                                         *
 * Maximum global history length and minimum global history length used    */

#ifndef MAXHIST
#define MAXHIST 131
#endif
#ifndef MINHIST
#define MINHIST 5
#endif



#define PERCEPTRON_HISTORY	62
#define NUM_PERCEPTRONS		2731
#define PERCEPTRON_BITS		8
#define MAX_WEIGHT		((1<<(PERCEPTRON_BITS-1))-1)
#define MIN_WEIGHT		(-(MAX_WEIGHT+1))
#define THETA			((int) (1.93 * PERCEPTRON_HISTORY + 14))
//#define THETA			((int) (1.93 * MAXHIST + 14))
/* size of buffer for keeping 'perceptron_state' for update */
#define NUM_UPDATE_ENTRIES	100



/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  END OF MACROS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


using namespace std;

typedef uint32_t address_t;
typedef bitset <MAXHIST> history_t;


/* This class defines cyclic shift register for folding long global history bits to smaller number of bits   * 
 * Original length of history bits for TAGE tables form a geometric series                                   *
 * Compressed length if used for indexing will then be LOGG number of bits                                   *
 * Compressed length if used for tag match, depends on number of tag bits for a given table                  */

class folded_history
{
public:
  unsigned comp;
  int CLENGTH;       // Compressed Length
  int OLENGTH;       // Original Length
  int OUTPOINT;      

  folded_history ()
  {
  }

  void init (int original_length, int compressed_length)
  {
    
    comp = 0;
    OLENGTH = original_length;
    CLENGTH = compressed_length;
    OUTPOINT = OLENGTH % CLENGTH;
    ASSERT (OLENGTH < MAXHIST);
  
  }


  /*Explain detail history folding*/

  void update (history_t h)
  {
    
    ASSERT ((comp >> CLENGTH) == 0);
    comp = (comp << 1) | h[0];
    comp ^= h[OLENGTH] << OUTPOINT;
    comp ^= (comp >> CLENGTH);
    comp &= (1 << CLENGTH) - 1;
  
  }
};


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  ALL THE PREDICTOR IS HERE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

class PREDICTOR
{
public:

  /* Bimodal table entry */
  class bentry
  {
  public:
    int8_t hyst;
    int8_t pred;

    /* Bimodal table consists of a 2 bit prediction counter (pred), initialized to 0   *
     * and a hysteresis bit which is shared among 4 entries, initialized to 1          *
     * pred (MSB of 2 bit counter) and hyst are either 0 or 1 always (1 bit counter)   */ 
    
    bentry ()
    {

      pred = 0;
      hyst = 1;

    }
  };

  
  /* Global TAGE table entry */
  class gentry       
  {
  public:
    int8_t ctr;
    uint16_t tag;
    int8_t ubit;

    /* TAGE table consists of a 3-bit prediction counter (ctr)                         *
     * variable length tag field, with a maximum of TBITS                              *
     * and a useful 2-bit entry (ubit) used for allocating new entry on update         *
     * All these fields are initialized to zero                                        */
    
    gentry ()
    {

      ctr = 0;       //~3 bit
      tag = 0;       //~9bit 
      ubit = 0;      //~2bit

    }
  };

  



 
 /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  PREDICTOR STRORAGE DATA - TAGE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  
  
  int PWIN;    /* 4 bit entry to determine newly allocated entries as valid or not */
  
  int TICK;    /* Utility variable for ageing useful counter for TAGE tables */
  
  int phist;   /* Use a path history for index computation as in the OGEHL predictor */
  
  
  history_t ghist;			/* CSR containing global history bits */	
  
  folded_history ch_i[NHIST];		/* Utility for computing TAGE indices */
  
  folded_history ch_t[2][NHIST];	/* Utility for computing TAGE tags */
  
  bentry *btable;			/* Bimodal TAGE table */
  
  gentry *gtable[NHIST];		/* Tagged TAGE tables */

  int m[NHIST];                         /* Used for storing the geometric history lengths */
  

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
    int STORAGESIZE = 0;
    ghist = 0;
    
    /* compute Geometric history lengths:                                                                      *
     * min * [ ( max / min ) ^ ( 1 / ( n-1 ) )  ] + 0.5                                                        *
     * This displaces intermediate values in a geometric series fashion with MINHIST and MAXHIST as bounds    */

    m[0] = MAXHIST - 1;
    m[NHIST - 1] = MINHIST;
    for (int i = 1; i < NHIST - 1; i++)
    {
      m[NHIST - 1 - i] = (int) (((double) MINHIST *
		  		pow ((double) (MAXHIST - 1) / (double) MINHIST,
		       			(double) (i) / (double) ((NHIST - 1)))) + 0.5);

    }

    //fprintf (stderr, "History Series:");
    
    /* Initialize ch_i with OLENGTH as corresponding history length and CLENGTH as LOGG */
    for (int i = NHIST - 1; i >= 0; i--)
    {

	//fprintf (stderr, "%d ", m[i]);
	ch_i[i].init (m[i], (LOGG));
	STORAGESIZE += (1 << LOGG) * (5 + TBITS - ((i + (NHIST & 1)) / 2));
    }
    //fprintf (stderr, "\n");
    STORAGESIZE += (1 << LOGB) + (1 << (LOGB - 2));
    fprintf (stderr, "NHIST= %d; MINHIST= %d; MAXHIST= %d; STORAGESIZE= %d bits\n",NHIST, MINHIST, MAXHIST - 1, STORAGESIZE);

    
    /* Initialize ch_t with OLENGTH same as ch_i and CLENGTH as below */
    for (int i = 0; i < NHIST; i++)
    {
	ch_t[0][i].init (ch_i[i].OLENGTH, TBITS - ((i + (NHIST & 1)) / 2));
	ch_t[1][i].init (ch_i[i].OLENGTH, TBITS - ((i + (NHIST & 1)) / 2) - 1);
    }

    btable = new bentry[1 << LOGB];          /* Bimodal table gets 2^LOGB entries */
    
    for (int i = 0; i < NHIST; i++)
      gtable[i] = new gentry[1 << (LOGG)];   /* Each of tagged tables gets 2^LOGG entries */
      
    /*******!!!! PERCEPTRON !!!!********/
    spec_global_history = 0;
    global_history = 0;
    perceptron_state_buf_ctr = 0;
    initialize_perceptron_predictor();


 }
/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  END OF CONSTRUCTOR >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/


  /* Indexes to the different tables are computed only once and store in GI*/
  int GI[NHIST];

  
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  F  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  
  /* Index function for the global tables:                     *
   * Includes path history as in the OGEHL predictor           *
   * Essentially mixes the path history                        */
  
  int F (int A, int size, int bank)
  {
    int A1, A2;

    A = A & ((1 << size) - 1);
    A1 = (A & ((1 << LOGG) - 1));
    A2 = (A >> LOGG);
    A2 = ((A2 << bank) & ((1 << LOGG) - 1)) + (A2 >> (LOGG - bank));
    A = A1 ^ A2;
    A = ((A << bank) & ((1 << LOGG) - 1)) + (A >> (LOGG - bank));
    return (A);
  }
 


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  gindex  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  /* Index function for the global tables:            *
   * Index = XOR ( PC, folded_hitory, path_history )  */
  
  int gindex (address_t pc, int bank)
  {
    int index;
    if (m[bank] >= 16)
      index = pc ^ (pc >> ((LOGG - (NHIST - bank - 1)))) ^ ch_i[bank].comp ^ F(phist, 16, bank);
    
    else 
      index = pc ^ (pc >> (LOGG - NHIST + bank + 1)) ^ ch_i[bank].comp ^ F(phist, m[bank], bank);

    
    return (index & ((1 << (LOGG)) - 1));

  }

 
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  gtag  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  /* Tag computation for global tables  * 
   * Tag = PC ^ CSR1 ^ (CSR2 << 1)      */
  
  uint16_t gtag (address_t pc, int bank)
  {

    int tag = pc ^ ch_t[0][bank].comp ^ (ch_t[1][bank].comp << 1);
    
    return (tag & ((1 << (TBITS - ((bank + (NHIST & 1)) / 2))) - 1));
    //return (tag & ((1 << (TBITS )) - 1));

  }
    
 
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  ctrupdate  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  /* Up-Down saturating counter update after prediction on actual branch outcome   *
   * nbits is number of bits for ctr (here 3)                                      *
   * 0,1,2 ... 7 is mapped to -4, -3, -2 ... 3                                   */

  void ctrupdate (int8_t & ctr, bool taken, int nbits)
  {
    if (taken)
      {
	if (ctr < ((1 << (nbits - 1)) - 1))    /* if ctr < 3 */
	  {
	    ctr++;
	  }
      }
    else
      {
	if (ctr > -(1 << (nbits - 1)))         /* if ctr > -4 */ 
	  {
	    ctr--;
	  }
      }
  }
  
  int altbank; /* Table corresponding to alternate prediction */
  
 

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  read_prediction  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  /* Prediction given by longest matching global history                                                 *
   * altpred contains the alternate prediction, second longest matching history, if none then bimodal    */

  bool read_prediction (address_t pc, int &bank, bool & altpred)
  {

    /* Initializing prediction and alternate prediction banks to NHIST */
    bank = NHIST;
    altbank = NHIST;
        

    /* Finding longest matching history */
    for (int i = 0; i < NHIST; i++)
    {
      if (gtable[i][GI[i]].tag == gtag (pc, i))
      {
        bank = i;
        break;
      }
    }
    

    /* Next longest matching history for altpred */
    for (int i = bank + 1; i < NHIST; i++)
    {
      if (gtable[i][GI[i]].tag == gtag (pc, i))
      {
        altbank = i;
        break;
      }
    }


    /* bank < NHIST => There was some tag match, else bimodal gives prediction and alternate prediction */
    if (bank < NHIST)
    {
      
      if (altbank < NHIST)
        altpred = (gtable[altbank][GI[altbank]].ctr >= 0);
      else
        altpred = getperceptron (pc);
        //altpred = getbim (pc);

      
      /* If prediction is given by newly allocated entry and if PWIN is >=0, then altpred is taken    *
       * Else if PWIN is -ve or not newly allocated entry prediction by longest match is taken        *
       * PWIN is a 4 bit counter to measure usefulness of altpred                                     */

      if ((PWIN < 0) || (abs (2 * gtable[bank][GI[bank]].ctr + 1) != 1) || (gtable[bank][GI[bank]].ubit != 0))
	return (gtable[bank][GI[bank]].ctr >= 0);
      else
        return (altpred);

    }
    else
    {
      //altpred = getbim (pc);
      altpred = getperceptron(pc);
      
      return altpred;
    }
    
  }//read_prediction


  bool pred_taken, alttaken;
  int bank;

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  get_prediction  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  /* PREDICTION by our predictor */
  
  bool get_prediction (const branch_record_c * br, const op_state_c * os)
  {
  
    if (br->is_conditional)
    {

	address_t pc = br->instruction_addr;     /* Address of the corresponding branch instruction */
 
        
	/* Indexes of global tagged tables and bimodal tables */
	for (int i = 0; i < NHIST; i++)
	  GI[i] = gindex (pc, i);
	
	BI = bindex (pc);

        
	/* bank contains the number of the matching table, NHIST if no match   *
         * pred_taken is the prediction                                        *
         * alttaken is the alternate prediction                                */
	
	pred_taken = read_prediction (pc, bank, alttaken);

    }
    
    return pred_taken;
  
  }//get_prediction
  

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  getperceptron  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  /* Prediction by the bimodal predictor */
  
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



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  MYRANDOM  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  /* Pseudo random number generator, generates alternate odd even values */
  
  int Seed;
  int MYRANDOM ()
  {
    Seed = ((1 << 2 * NHIST) + 1) * Seed + 0xf3f531;
    Seed = (Seed & ((1 << (2 * (NHIST))) - 1));
    return (Seed);
  };
 

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  update_predictor  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  /* Update the predictor after the branch outcome is known comparing with the prediction by our predictor */

  void update_predictor (const branch_record_c * br, const op_state_c * os, bool taken)
  {

    int NRAND = MYRANDOM ();

    /* We predict only for conditional branches, i.e direct branches */
    if (br->is_conditional)
    {
	
	address_t pc = br->instruction_addr;

	/* We allocate a new entry in strictly one of the tagged tables if prediction by predictor is wrong and    *
	 * provider component is not the one with longest history                                                  */	
		bool ALLOC = ((pred_taken != taken) & (bank > 0));


	/* On a TAG hit, we find the usefulness of alternate prediction if it differs from predictor component     *
	 * Monitor the usefulness through a 4-bit counter PWIN                                                     */
	
	if (bank < NHIST)
	{
	    
	    /* True if not a new entry and predicted taken */ 
	    bool loctaken = (gtable[bank][GI[bank]].ctr >= 0);  
	    
	    /* True if the prediction is by newly allocated entry, ctr value for a new entry can be -1 or 1 and ubit is 0 */
	    bool PseudoNewAlloc = (abs (2 * gtable[bank][GI[bank]].ctr + 1) == 1) && (gtable[bank][GI[bank]].ubit == 0);


	    if (PseudoNewAlloc)
	    {

		/* If the provider component was delivering the correct prediction; no need to allocate a new entry    * 
		 * even if the overall prediction was false                                                             */
		
		if (loctaken == taken)
		  ALLOC = false;


		if (loctaken != alttaken)   /* altpred differs from pred */
		{
		    if (alttaken == taken)  /* altpred correctly predicted */ 
		    {

			if (PWIN < 7)
			  PWIN++;
		    }

		    else if (PWIN > -8)
		      PWIN--;
		}
	    }
	}//if (bank < NHIST)

	
	/* Allocate a new entry */
	if (ALLOC)
	{

            /* Is there some "unuseful" entry to evacuate and allocate new???          *
	     * An entry is deemed "unuseful" if corresponding ubit is zero             *
	     * (i.e) it was newly allocated and never provided a correct prediction    */

	    int8_t min = 3;
	    for (int i = 0; i < bank; i++)
	    {
		
		if (gtable[i][GI[i]].ubit < min)
		  min = gtable[i][GI[i]].ubit;

	    }

	    
            /* No UNUSEFUL entry to allocate: Age all possible targets, but do not allocate */
	    if (min > 0)
	    {

		/* Ageing is done for tables with tables having larger history lengths than the predictor component */
		for (int i = bank - 1; i >= 0; i--)
		{
		    gtable[i][GI[i]].ubit--;

		}

	    }
            /* YES: Allocate one entry, but apply some randomness   *
             * bank I is twice more probable than bank I-1          */
	    else
	    {

		/* For every alternate iteration, NRAND is odd and bank I (smaller history) is chosen over I-1 if ubit is 0 */
		int Y = NRAND & ((1 << (bank - 1)) - 1);
		int X = bank - 1;
		while ((Y & 1) != 0) //Loops as long as bits in Y are 1
		{
		    
		    X--;
		    Y >>= 1;

		}

                
		for (int i = X; i >= 0; i--)
		{
		    int T = i;
		    if ((gtable[T][GI[T]].ubit == min)) 
		    {

			gtable[T][GI[T]].tag = gtag (pc, T);       //Corresponding TAG that was not matched
			gtable[T][GI[T]].ctr = (taken) ? 0 : -1;   //If taken ctr is 0 (weak-taken), else -1 (weak-not taken) 
			gtable[T][GI[T]].ubit = 0;                 //usefulness bit is set to 0
			break;
		    }
		}

	    }

	}//if(ALLOC)



/* Periodic reset of ubit (2-bit counter) every 256k branches: reset is not complete but bit by bit          *
 * First MSB is reset, next time LSB is reset and this iterates over every 256k branch instructions once     *
 * For once in two 256k, 18th bit becomes 1. And x is set 2 which resets LSB,                                *  
 * Else x is 1 and MSB is reset across all entries @ all tables                                              */  
	
	TICK++;
	if ((TICK & ((1 << 18) - 1)) == 0)     // 1 << 18 is 256k
	{
	    
	    int X = (TICK >> 18) & 1;
	    
	    if ((X & 1) == 0)
	      X = 2;
	    
	    for (int i = 0; i < NHIST; i++)
	      for (int j = 0; j < (1 << LOGG); j++)
		gtable[i][j].ubit = gtable[i][j].ubit & X;
	  
	}

	
	/*  Update the counter that provided the prediction, and only that counter [Tag tables or Bimodal]    *
	    Atmost one update of ctr (3-bit counter)                                                          *
	    Decrement of counters @ other tables not effected as Tag match might not be ther	              */
	
	if (bank < NHIST) 
	    ctrupdate (gtable[bank][GI[bank]].ctr, taken, CBITS);
	else
            perceptronupdate (pc, taken);
	
	
	/*  Useful counter of provider component updated only when alt.prediction differs final prediction    *
	    ubit (2-bit counter) incremented on taken branch and decremented on not-taken                     */
	
	if ((pred_taken != alttaken))
	{
	    ASSERT (bank < NHIST);

	    if (pred_taken == taken)
	    {

		if (gtable[bank][GI[bank]].ubit < 3)
		  gtable[bank][GI[bank]].ubit++;

	    }
	    else
	    {
		
		if (gtable[bank][GI[bank]].ubit > 0)
		  gtable[bank][GI[bank]].ubit--;
	      
	    }

	}

    }


/* Update global history and path history CSRs                                   *
 * Use also histoy bits of unconditional branches as for OGEHL predictors        *
 * Path history takes one bit (LSB) per branch-instruction pc                    *  
 * Path history is a 16 bit CSR                                                  */  

    ghist = (ghist << 1);
    if ((!br->is_conditional) | (taken))
      ghist |= (history_t) 1;

    phist = (phist << 1) + (br->instruction_addr & 1);
    phist = (phist & ((1 << 16) - 1));
    
    
    /* After updating history, history folding @ geometric history length tables happen */
    
    for (int i = 0; i < NHIST; i++)
    {
	ch_i[i].update (ghist);
	ch_t[0][i].update (ghist);
	ch_t[1][i].update (ghist);
    }

 }//update_predictor

};
#endif // PREDICTOR_H_SEEN

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  END ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
