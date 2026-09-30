#include "m_pd.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

// second-order markov chain generator
// needs a sequence of at least 3 values as data (set with message e.g. "[setseq 69 70 71 69 69 71( ")
// optionally set the starting state with a float (set with message e.g. "[setstart 71( ")
// bang to generate next value

static t_class *markov_class;

typedef struct _entry
{
	t_float values[3]; // prevprev, prev, curr
	t_float prob;			 // probability
} t_entry;

typedef struct _matrix
{
	t_entry *m; // entries
	int size;		//
} t_matrix;

typedef struct _markov
{
	t_object obj;					 // the Pd object
	t_float *seq;					 // array of floats representing the sequence to be turned into a probability matrix
	int seq_size;					 // size of the above array
	t_float prevprev;			 // value before the value before the previous
	t_float prev;					 // previous value
	t_float curr;					 // current value
	t_matrix *prob_matrix; // probability matrix
	t_outlet *out;
} t_markov;

// patch-level methods
void *markov_new(t_symbol *n);
void markov_free(t_markov *x);
void markov_bang(t_markov *x);
void markov_setseq(t_markov *x, t_symbol *s, int argc, t_atom *argv);
void markov_setstart(t_markov *x, t_floatarg f);

// helper methods
void prob_matrix_set(t_markov *x, t_matrix *mat, t_float pp, t_float p, t_float c, t_float prob);
float prob_matrix_get(t_markov *x, t_matrix *mat, t_float pp, t_float p, t_float c);
void print_matrix(t_markov *x);
int entry_cmp(const void *a, const void *b);

void markov_setup(void)
{
	t_class *c;
	markov_class = class_new(gensym("markov"), (t_newmethod)markov_new, (t_method)markov_free,
													 sizeof(t_markov), 0, A_NULL, 0);
	c = markov_class;

	class_addbang(c, (t_method)markov_bang);
	class_addmethod(c, (t_method)markov_setstart, gensym("setstart"), A_FLOAT, 0);
	class_addmethod(c, (t_method)markov_setseq, gensym("setseq"), A_GIMME, 0);
}

void *markov_new(t_symbol *myname)
{
	t_markov *x = (t_markov *)pd_new(markov_class);
	x->seq = NULL;
	x->seq_size = 0;
	x->prevprev = -1;
	x->prev = -1;
	x->curr = -1;
	x->prob_matrix = (t_matrix *)malloc(sizeof(t_matrix));
	x->prob_matrix->m = NULL;
	x->prob_matrix->size = 0;
	x->out = outlet_new(&x->obj, &s_float);

	srand((unsigned)time(NULL));
	return (void *)x;
}

void markov_free(t_markov *x)
{
	if (x->seq)
	{
		free(x->seq);
		x->seq = NULL;
	}
	if (x->prob_matrix)
	{
		free(x->prob_matrix->m);
		x->prob_matrix->m = NULL;
		free(x->prob_matrix);
		x->prob_matrix = NULL;
	}
}

// called when setseq <float array> sent to markov
// sets the sequence, generates the probability matrix, then prints it onto the console
// outputs an error if the given sequence does not have at least 3 values, or if any of the values are not floats >= 0
// treats the sequence as circular, e.g. "[setseq 61 62 65 63(" sets 63 61 as a valid pair, with 62 as the possible next value.
void markov_setseq(t_markov *x, t_symbol *s, int argc, t_atom *argv)
{
	// free memory to set a new sequence and probability matrix
	if (x->seq != NULL)
	{
		free(x->seq);
		x->seq = NULL;
	}
	if (x->prob_matrix != NULL)
	{
		free(x->prob_matrix->m);
		x->prob_matrix->m = NULL;
		free(x->prob_matrix);
		x->prob_matrix = NULL;
	}

	// check to make sure there are at least 3 values given
	if (argc > 2)
	{
		x->seq = (t_float *)malloc(argc * sizeof(t_float));
		if (!x->seq)
		{
			pd_error(x, "markov: setseq -- out of memory");
			x->seq_size = 0;
			return;
		}
		x->seq_size = argc;

		// ensure each element is a float >= 0
		for (int i = 0; i < argc; i++)
			if (argv[i].a_type == A_FLOAT && atom_getfloat(argv + i) >= 0.0f)
				x->seq[i] = atom_getfloat(argv + i);
			else if (argv[i].a_type != A_FLOAT)
			{
				x->seq[i] = 0.0f;
				pd_error(x, "markov: setseq -- element %d is not a float", i);
			}
			else
			{
				x->seq[i] = 0.0f;
				pd_error(x, "markov: setseq -- element %d is not within the accepted range", i);
			}

		x->prob_matrix = (t_matrix *)malloc(sizeof(t_matrix));
		x->prob_matrix->m = NULL;
		x->prob_matrix->size = 0;

		t_matrix *count_matrix = (t_matrix *)malloc(sizeof(t_matrix));
		count_matrix->m = NULL;
		count_matrix->size = 0;

		// count occurences of each possible next value of a given pair
		// separately count the total number of times that pair occurs
		// finally, normalize by dividing each possible next value of a given pair
		// by the number of time that pair occurs
		t_float pp = x->seq[argc - 2];
		t_float p = x->seq[argc - 1];
		t_float c = x->seq[0];
		t_float get_total = prob_matrix_get(x, count_matrix, pp, p, -1);
		t_float get = prob_matrix_get(x, x->prob_matrix, pp, p, c);
		prob_matrix_set(x, count_matrix, pp, p, -1, get_total + 1);
		prob_matrix_set(x, x->prob_matrix, pp, p, c, get + 1);

		pp = x->seq[argc - 1];
		p = x->seq[0];
		c = x->seq[1];
		get_total = prob_matrix_get(x, count_matrix, pp, p, -1);
		get = prob_matrix_get(x, x->prob_matrix, pp, p, c);
		prob_matrix_set(x, count_matrix, pp, p, -1, get_total + 1);
		prob_matrix_set(x, x->prob_matrix, pp, p, c, get + 1);

		for (int i = 2; i < argc; i++)
		{
			pp = x->seq[i - 2];
			p = x->seq[i - 1];
			c = x->seq[i];
			get_total = prob_matrix_get(x, count_matrix, pp, p, -1);
			get = prob_matrix_get(x, x->prob_matrix, pp, p, c);
			prob_matrix_set(x, count_matrix, pp, p, -1, get_total + 1);
			prob_matrix_set(x, x->prob_matrix, pp, p, c, get + 1);
		}

		// normalize
		for (int i = 0; i < x->prob_matrix->size; i++)
		{
			t_float epp = x->prob_matrix->m[i].values[0];
			t_float ep = x->prob_matrix->m[i].values[1];
			t_float total = prob_matrix_get(x, count_matrix, epp, ep, -1);
			if (total > 0)
				x->prob_matrix->m[i].prob /= total;
		}

		free(count_matrix->m);
		count_matrix->m = NULL;
		free(count_matrix);
		count_matrix = NULL;

		// print
		print_matrix(x);
	}
	else
	{
		pd_error(x, "markov: setseq -- sequence must have 3 or more elements");
		x->seq_size = 0;
	}
}

// called when setstart <float> sent to markov
// sets the starting state of the markov chain
// outputs an error if the sequence is not already set, or if the given float is not within the sequence
// otherwise, always finds the first instance of the float within the sequence
// the two indicies in the sequence before first instance of the given float are used to set the initial state
void markov_setstart(t_markov *x, t_floatarg f)
{
	if (!x->seq || x->seq_size <= 0)
	{
		pd_error(x, "markov: setstart -- no sequence defined, send 'setseq' symbol with sequence (e.g. '[setseq 69 70 71 69 69 71(' ) first");
		return;
	}

	// look for given float
	int found = 0;
	for (int i = 0; i < x->seq_size; i++)
		if (x->seq[i] == f)
		{
			found = 1;
			break;
		}

	// if found, set to starting state to first instance of given float and the two indices before that in the sequence
	if (found)
	{
		x->curr = f;
		int idx = 0;
		for (int i = 0; i < x->seq_size; i++)
			if (x->curr == x->seq[i])
			{
				idx = i;
				break;
			}
		if (idx >= 2)
		{
			x->prevprev = x->seq[idx - 2];
			x->prev = x->seq[idx - 1];
		}
		else if (idx == 1)
		{
			x->prevprev = x->seq[x->seq_size - 1];
			x->prev = x->seq[0];
		}
		else
		{
			x->prevprev = x->seq[x->seq_size - 2];
			x->prev = x->seq[x->seq_size - 1];
		}
		post("markov: setstart -- start value set to %i, with previous and doubly-previous values established as %i and %i, respectively.", (int)x->curr, (int)x->prev, (int)x->prevprev);
	}
	else
		pd_error(x, "markov: setstart -- value %i is not in the currently loaded sequence", (int)f);
}

// called when <bang> sent to markov
// outputs the next value based on markov chain probabilities and updates the state
// outputs an error if the starting sequence was not defined
void markov_bang(t_markov *x)
{
	if (x->prob_matrix && x->seq && x->seq_size > 0)
	{
		// if no starting state was given, choose the 1st, 2nd, and 3rd floats to define the starting state
		if (x->curr == -1)
		{
			x->prevprev = x->seq[0];
			x->prev = x->seq[1];
			x->curr = x->seq[2];
		}

		outlet_float(x->obj.ob_outlet, x->curr);

		// use a random selection with weights based on the probabilities from the probability matrix
		// to select the next value, then update the state
		x->prevprev = x->prev;
		x->prev = x->curr;
		t_float next = x->curr;
		t_float cumulative_prob = 0.0f;
		t_float rand_float = (t_float)rand() / (t_float)RAND_MAX;

		for (int i = 0; i < x->prob_matrix->size; i++)
		{
			if (x->prob_matrix->m[i].values[0] == x->prevprev &&
					x->prob_matrix->m[i].values[1] == x->prev)
			{
				cumulative_prob += x->prob_matrix->m[i].prob;
				if (cumulative_prob >= rand_float)
				{
					next = x->prob_matrix->m[i].values[2];
					break;
				}
			}
		}

		x->curr = next;
	}
	else
		post("markov: bang -- no sequence defined, send 'setseq' symbol with sequence (e.g. '[setseq 69 70 71 69 69 71(' ) first");
}

// helper function that sets values in the probability matrix
// or, if the given key does not already exist within the matrix, adds that entry
void prob_matrix_set(t_markov *x, t_matrix *mat, t_float pp, t_float p, t_float c, t_float prob)
{
	// search for existing key to update
	for (int i = 0; i < mat->size; i++)
		if (mat->m[i].values[0] == pp &&
				mat->m[i].values[1] == p &&
				mat->m[i].values[2] == c)
		{
			mat->m[i].prob = prob;
			return;
		}

	// expand the array using realloc if key not found
	int new_size = mat->size + 1;
	t_entry *temp = (t_entry *)realloc(mat->m, new_size * sizeof(t_entry));

	// check if out of memory
	if (temp == NULL)
	{
		pd_error(x, "markov: prob_matrix_set -- map resize failed, out of memory");
		return;
	}

	// assign new memory and add entry
	mat->m = temp;
	mat->m[mat->size].values[0] = pp;
	mat->m[mat->size].values[1] = p;
	mat->m[mat->size].values[2] = c;
	mat->m[mat->size].prob = prob;
	mat->size = new_size;
}

// helper function to get probability from a probability matrix
// initializes an entry with value 0 and given values as key if key does not exist within matrix
float prob_matrix_get(t_markov *x, t_matrix *mat, t_float pp, t_float p, t_float c)
{
	for (int i = 0; i < mat->size; i++)
	{
		if (mat->m[i].values[0] == pp &&
				mat->m[i].values[1] == p &&
				mat->m[i].values[2] == c)
		{
			return mat->m[i].prob;
		}
	}

	prob_matrix_set(x, mat, pp, p, c, 0.0f);
	return 0.0f;
}

// helper function to print the probability matrix
void print_matrix(t_markov *x)
{
	t_entry *sorted = (t_entry *)malloc(x->prob_matrix->size * sizeof(t_entry));
	if (!sorted)
	{
		pd_error(x, "markov: print_matrix -- out of memory");
		return;
	}

	// create a copy of the probability matrix, then sort the copy
	memcpy(sorted, x->prob_matrix->m, x->prob_matrix->size * sizeof(t_entry));
	qsort(sorted, x->prob_matrix->size, sizeof(t_entry), entry_cmp);

	post("markov: probability matrix with %d entries", x->prob_matrix->size);
	post("  %-22s %-10s %-10s", "-----prev/curr-----", "next", "prob");
	post("  %-22s %-10s %-10s", "-------------------", "----", "----");

	t_float cur_pp = sorted[0].values[0];
	t_float cur_p = sorted[0].values[1];

	for (int i = 0; i < x->prob_matrix->size; i++)
	{
		t_float pp = sorted[i].values[0];
		t_float p = sorted[i].values[1];
		t_float c = sorted[i].values[2];
		t_float prob = sorted[i].prob;

		if (pp != cur_pp || p != cur_p)
		{
			post("");
			cur_pp = pp;
			cur_p = p;
		}

		post("  %-10.2f %-10.2f %-10.2f %.4f", pp, p, c, prob);
	}
	post("");

	free(sorted);
}

// comparator function to help sort the entries in the probability matrix for more interpretable printing
int entry_cmp(const void *a, const void *b)
{
	const t_entry *ea = (const t_entry *)a;
	const t_entry *eb = (const t_entry *)b;

	if (ea->values[0] != eb->values[0])
		return (ea->values[0] > eb->values[0]) - (ea->values[0] < eb->values[0]);
	if (ea->values[1] != eb->values[1])
		return (ea->values[1] > eb->values[1]) - (ea->values[1] < eb->values[1]);
	return (ea->values[2] > eb->values[2]) - (ea->values[2] < eb->values[2]);
}