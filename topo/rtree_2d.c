/****************************************************************************
* RTree.C
*
* MODULE:       R-Tree library 
*              
* AUTHOR(S):    Antonin Guttman - original code
*               Daniel Green (green@superliminal.com) - major clean-up
*                               and implementation of bounding spheres
*               
* PURPOSE:      Multi Dimensional Index
*
* COPYRIGHT:    (C) 2001 by the GRASS Development Team
*
*               This program is free software under the GNU General Public
*               License (>=v2). Read the file COPYING that comes with GRASS
*               for details.
*
* LAST MODIFY:         ZhangLiang (cheungmine@gmail.com) - 2007-11
*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <math.h>

#include "rtree_2d.h"

#define		METHODS		1

/* variables for finding a partition */
typedef struct _RTREEPARTITION
{
	int			partition[MAXCARD+1];
	int			total;
	int			minfill;
	int			taken[MAXCARD+1];
	int			count[2];
	RTREEMBR	cover[2];
	REALTYPE	area[2];
} RTREEPARTITION;

typedef struct _RTREENODE
{
	int	count;
	int	level; /* 0 is leaf, others positive */
	RTREEBRANCH  branch[MAXCARD];
}RTREENODE;

typedef struct _RTREELISTNODE
{
     struct _RTREELISTNODE	*next;
     RTREENODE		*node;
}RTREELISTNODE;

typedef struct _RTREEROOT
{
	RTREENODE*		root_node;
	
	RTREEBRANCH		BranchBuf[MAXCARD+1];
	int				BranchCount;
	RTREEMBR		CoverSplit;
	REALTYPE		CoverSplitArea;
	RTREEPARTITION	Partitions[METHODS];

	PFNRTreeSearchCallback  pfnCallback;
} RTREEROOT;


#define BIG_NUM (FLT_MAX/4.0)

#define INVALID_RECT(x) ((x)->bound[0] > (x)->bound[DIMS_NUMB])

#ifndef MIN
	#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
	#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

int NODECARD = MAXCARD;
int LEAFCARD = MAXCARD;

/* balance criteria for node splitting */
/* NOTE: can be changed if needed. */
#define MINNODEFILL (NODECARD / 2)
#define MINLEAFFILL (LEAFCARD / 2)

#define MAXKIDS(n) ((n)->level > 0 ? NODECARD : LEAFCARD)
#define MINFILL(n) ((n)->level > 0 ? MINNODEFILL : MINLEAFFILL)

/**
 * Methods below not used at all!
 */
static int set_max(int *which, int new_max){
	if(2 > new_max || new_max > MAXCARD)
		return 0;
	*which = new_max;
	return 1;
}
static int RTreeSetNodeMax(int new_max) { return set_max(&NODECARD, new_max); }
static int RTreeSetLeafMax(int new_max) { return set_max(&LEAFCARD, new_max); }
static int RTreeGetNodeMax(void) { return NODECARD; }
static int RTreeGetLeafMax(void) { return LEAFCARD; }


/**
 * Load branch buffer with branches from full node plus the extra branch.
 */
static void _RTreeGetBranches(HRTREEROOT root, RTREENODE *node, RTREEBRANCH *br)
{
	int i;

	assert(node && br);
	
	/* load the branch buffer */
	for (i=0; i<MAXKIDS(node); i++)
	{
		assert(node->branch[i].child); /* n should have every entry full */
		root->BranchBuf[i] = node->branch[i];
	}

	root->BranchBuf[MAXKIDS(node)] = *br;
	root->BranchCount = MAXKIDS(node) + 1;

	/* calculate mbr containing all in the set */
	root->CoverSplit = root->BranchBuf[0].mbr;

	for (i=1; i<MAXKIDS(node)+1; i++)
	{
		root->CoverSplit = RTreeCombineRect(&root->CoverSplit, &root->BranchBuf[i].mbr);
	}

	root->CoverSplitArea = RTreeRectSphericalVolume(&root->CoverSplit);
	RTreeInitNode(node);
}



/**
 * Put a branch in one of the groups.
 */
static void _RTreeClassify(HRTREEROOT root, int i, int group, RTREEPARTITION *p)
{
	assert(p);
	assert(!p->taken[i]);

	p->partition[i] = group;
	p->taken[i] = TRUE;

	if (p->count[group] == 0)
		p->cover[group] = root->BranchBuf[i].mbr;
	else
		p->cover[group] = RTreeCombineRect(&root->BranchBuf[i].mbr, &p->cover[group]);
	
	p->area[group] = RTreeRectSphericalVolume(&p->cover[group]);
	p->count[group]++;
}

/**
 * Pick two rects from set to be the first elements of the two groups.
 * Pick the two that waste the most area if covered by a single rectangle.
 */
static void _RTreePickSeeds(HRTREEROOT root, RTREEPARTITION *p)
{
	int i, j, seed0=0, seed1=0;
	REALTYPE worst, waste, area[MAXCARD+1];

	for (i=0; i<p->total; i++)
		area[i] = RTreeRectSphericalVolume(&root->BranchBuf[i].mbr);
	
	worst = -root->CoverSplitArea - 1;
	
	for (i=0; i<p->total-1; i++)
	{
		for (j=i+1; j<p->total; j++)
		{
			RTREEMBR one_rect;
			one_rect = RTreeCombineRect(&root->BranchBuf[i].mbr, &root->BranchBuf[j].mbr);
			waste = RTreeRectSphericalVolume(&one_rect) - area[i] - area[j];
			if (waste > worst)
			{
				worst = waste;
				seed0 = i;
				seed1 = j;
			}
		}
	}
	_RTreeClassify(root, seed0, 0, p);
	_RTreeClassify(root, seed1, 1, p);
}


/**
 * Copy branches from the buffer into two nodes according to the partition.
 */
static void _RTreeLoadNodes(HRTREEROOT root, RTREENODE *n, RTREENODE *q, RTREEPARTITION *p)
{
	int i;
	assert(n && q && p);

	for (i=0; i<p->total; i++)
	{
		assert(p->partition[i] == 0 || p->partition[i] == 1);
		if (p->partition[i] == 0)
			RTreeAddBranch(root, &root->BranchBuf[i], n, NULL);
		else if (p->partition[i] == 1)
			RTreeAddBranch(root, &root->BranchBuf[i], q, NULL);
	}
}

/**
 * Initialize a RTREEPARTITION structure.
 */
static void _RTreeInitPart( RTREEPARTITION *p, int maxrects, int minfill)
{
	int i;
	assert(p);

	p->count[0] = p->count[1] = 0;
	p->cover[0] = p->cover[1] = RTreeNullRect();
	p->area[0] = p->area[1] = (REALTYPE)0;
	p->total = maxrects;
	p->minfill = minfill;
	for (i=0; i<maxrects; i++)
	{
		p->taken[i] = FALSE;
		p->partition[i] = -1;
	}
}


/**
 * Print out data for a partition from RTREEPARTITION struct.
 */
static void _RTreePrintPart(HRTREEROOT root, RTREEPARTITION *p)
{
	int i;
	assert(p);
	
	fprintf (stdout, "\npartition:\n");
	for (i=0; i<p->total; i++)
	{
		fprintf (stdout, "%3d\t", i);
	}
	fprintf (stdout, "\n");
	for (i=0; i<p->total; i++)
	{
		if (p->taken[i])
			fprintf (stdout, "  t\t");
		else
			fprintf (stdout, "\t");
	}
	fprintf (stdout, "\n");
	for (i=0; i<p->total; i++)
	{
		fprintf (stdout, "%3d\t", p->partition[i]);
	}
	fprintf (stdout, "\n");

	fprintf (stdout, "count[0] = %d  area = %f\n", p->count[0], p->area[0]);
	fprintf (stdout, "count[1] = %d  area = %f\n", p->count[1], p->area[1]);
	if (p->area[0] + p->area[1] > 0)
	{
		fprintf (stdout, "total area = %f  effectiveness = %3.2f\n", 
			p->area[0] + p->area[1], (float)root->CoverSplitArea / (p->area[0] + p->area[1]));
	}
	fprintf (stdout, "cover[0]:\n");
	RTreePrintRect(&p->cover[0], 0);

	fprintf (stdout, "cover[1]:\n");
	RTreePrintRect(&p->cover[1], 0);
}


/**
 * Method #0 for choosing a partition:
 * As the seeds for the two groups, pick the two rects that would waste the
 * most area if covered by a single rectangle, i.e. evidently the worst pair
 * to have in the same group.
 * Of the remaining, one at a time is chosen to be put in one of the two groups.
 * The one chosen is the one with the greatest difference in area expansion
 * depending on which group - the mbr most strongly attracted to one group
 * and repelled from the other.
 * If one group gets too full (more would force other group to violate min
 * fill requirement) then other group gets the rest.
 * These last are the ones that can go in either group most easily.
 */
static void _RTreeMethodZero(HRTREEROOT root, RTREEPARTITION *p, int minfill )
{
	int i;
	REALTYPE biggestDiff;
	int group, chosen=0, betterGroup=0;
	assert(p);

	_RTreeInitPart(p, root->BranchCount, minfill);
	_RTreePickSeeds(root, p);

	while (p->count[0] + p->count[1] < p->total && 
		p->count[0] < p->total - p->minfill && 
		p->count[1] < p->total - p->minfill)
	{
		biggestDiff = (REALTYPE)-1.;
		for (i=0; i<p->total; i++)
		{
			if (!p->taken[i])
			{
				RTREEMBR *r, rect_0, rect_1;
				REALTYPE growth0, growth1, diff;
				
				r = &root->BranchBuf[i].mbr;
				rect_0 = RTreeCombineRect(r, &p->cover[0]);
				rect_1 = RTreeCombineRect(r, &p->cover[1]);
				growth0 = RTreeRectSphericalVolume(&rect_0) - p->area[0]; 
				growth1 = RTreeRectSphericalVolume(&rect_1) - p->area[1];
				diff = growth1 - growth0;
				if (diff >= 0)
					group = 0;
				else
				{
					group = 1;
					diff = -diff;
				}
				if (diff > biggestDiff)
				{
					biggestDiff = diff;
					chosen = i;
					betterGroup = group;
				}
				else if (diff==biggestDiff && p->count[group]<p->count[betterGroup])
				{
					chosen = i;
					betterGroup = group;
				}
			}
		}
		_RTreeClassify(root, chosen, betterGroup, p);
	}
	
	/* if one group too full, put remaining rects in the other */
	if (p->count[0] + p->count[1] < p->total)
	{
		if (p->count[0] >= p->total - p->minfill)
			group = 1;
		else
			group = 0;
		
		for (i=0; i<p->total; i++)
		{
			if (!p->taken[i])
				_RTreeClassify(root, i, group, p);
		}
	}
	
	assert(p->count[0] + p->count[1] == p->total);
	assert(p->count[0] >= p->minfill && p->count[1] >= p->minfill);
}


/**
 * Initialize one branch cell in a node. 
 */
static void _RTreeInitBranch( RTREEBRANCH *br )
{
	RTreeInitRect(&(br->mbr));
	br->child = NULL;
}


static void _RTreePrintBranch( RTREEBRANCH *br, int depth )
{
	RTreePrintRect(&(br->mbr), depth);
	RTreePrintNode(br->child, depth);
}


/**
 * Inserts a new data rectangle into the index structure.
 * Recursively descends tree, propagates splits back up.
 * Returns 0 if node was not split.  Old node updated.
 * If node was split, returns 1 and sets the pointer pointed to by
 * new_node to point to the new node.  Old node updated to become one of two.
 * The level argument specifies the number of steps up from the leaf
 * level to insert; e.g. a data rectangle goes in at level = 0.
 */
static int _RTreeInsertRect(HRTREEROOT root, RTREEMBR *mbr, void* tid,  RTREENODE *node, RTREENODE **new_node, int level)
{
	int i;
	RTREEBRANCH b;
	RTREENODE *n2;

	assert(mbr && node && new_node);
	assert(level >= 0 && level <= node->level);

	/* Still above level for insertion, go down tree recursively */
	if (node->level > level)
	{
		i = RTreePickBranch(mbr, node);
		if (!_RTreeInsertRect(root, mbr, tid, node->branch[i].child, &n2, level))
		{
			/* child was not split */
			node->branch[i].mbr = RTreeCombineRect(mbr, &(node->branch[i].mbr));
			return 0;
		}
		
		/* child was split */
		node->branch[i].mbr = RTreeNodeCover(node->branch[i].child);
		b.child = n2;
		b.mbr = RTreeNodeCover(n2);

		return RTreeAddBranch(root, &b, node, new_node);
	}	
	else if (node->level == level)	/* Have reached level for insertion. Add mbr, split if necessary */
	{
		b.mbr = *mbr;

#pragma warning(push)	/* C4312 */
#pragma warning( disable : 4312 )
		b.child = ( RTREENODE *) tid;
#pragma warning(pop)

		/* child field of leaves contains tid of data record */
		return RTreeAddBranch(root, &b, node, new_node);
	}
	
	/* Not supposed to happen */
	assert (FALSE);
	return 0;
}


/**
 * Allocate space for a node in the list used in DeletRect to
 * store Nodes that are too empty.
 */
static RTREELISTNODE * _RTreeNewListNode(void)
{
	return (RTREELISTNODE *) malloc(sizeof(RTREELISTNODE));
}

static void _RTreeFreeListNode(RTREELISTNODE *p)
{
	free(p);
}

/**
 * Add a node to the reinsertion list.  All its branches will later
 * be reinserted into the index structure.
 */
static void _RTreeReInsert(RTREENODE *node, RTREELISTNODE **ne)
{
	RTREELISTNODE *ln = _RTreeNewListNode();
	ln->node = node;
	ln->next = *ne;
	*ne = ln;
}

/**
 * Delete a rectangle from non-root part of an index structure.
 * Called by RTreeDeleteRect.  Descends tree recursively,
 * merges branches on the way back up.
 * Returns 1 if record not found, 0 if success.
 */
static int _RTreeDeleteRect( RTREEMBR *mbr, void* tid, RTREENODE *node, RTREELISTNODE **ee)
{
	int i;

	assert(mbr && node && ee);
	assert(node->level >= 0);

	if (node->level > 0)  /* not a leaf node */
	{
		for (i = 0; i < NODECARD; i++)
		{
			if (node->branch[i].child && RTreeOverlap( mbr, &(node->branch[i].mbr )))
			{
				if (!_RTreeDeleteRect( mbr, tid, node->branch[i].child, ee ))
				{
					if (node->branch[i].child->count >= MINNODEFILL) 
						node->branch[i].mbr = RTreeNodeCover(	node->branch[i].child );
					else{	/* not enough entries in child, eliminate child node */
						_RTreeReInsert(node->branch[i].child, ee);
						RTreeDisconnectBranch(node, i);
					}
					return 0;
				}
			}
		}
		return 1;
	}

#pragma warning(push)	/* C4312 */
#pragma warning( disable : 4312 )

	/* a leaf node */
	for (i = 0; i < LEAFCARD; i++)
	{
		if ( node->branch[i].child && node->branch[i].child == (RTREENODE *) tid )
		{
			RTreeDisconnectBranch( node, i );
			return 0;
		}

	}
#pragma warning(pop)

	return 1;
}

static void _RTreeTabIn(int depth)
{
	int i;
	for(i=0; i<depth; i++)
		putchar('\t');
}


/**
 * Search in an index tree or subtree for all data rectangles that overlap the argument rectangle.
 * Return the number of qualifying data rects.
 */
static int _RTreeSearchRect( RTREENODE *node, RTREEMBR *mbr, PFNRTreeSearchCallback pfnSHCB, void* pfnParam)
{
	/* Fix not yet tested. */
	int hitCount = 0;
	int i;

	assert(node && mbr);
	assert(node->level >= 0);
	
	if (node->level > 0) /* this is an internal node in the tree */
	{
		for (i=0; i<NODECARD; i++){
			if (node->branch[i].child && RTreeOverlap(mbr, &node->branch[i].mbr))
				hitCount += _RTreeSearchRect(node->branch[i].child, mbr, pfnSHCB, pfnParam);
		}
	}
	else /* this is a leaf node */
	{
#pragma warning(push)	/* C4311 */
#pragma warning( disable : 4311 )
		for (i=0; i<LEAFCARD; i++)
		{
			if (node->branch[i].child && RTreeOverlap(mbr, &node->branch[i].mbr))
			{
				hitCount++;

				/* call the user-provided callback and return if callback wants to terminate search early */
				if(pfnSHCB && ! pfnSHCB((void*)node->branch[i].child, pfnParam) )
					return hitCount; 

			}
		}
#pragma warning(pop)
	}
	return hitCount;
}

/*=============================================================================
								Public functions:
 =============================================================================*/
/**
 * Initialize a rectangle to have all 0 coordinates.
 */
 void RTreeInitRect( RTREEMBR *mbr)
{
	int i;
	for (i=0; i<SIDES_NUMB; i++)
		mbr->bound[i] = (REALTYPE) 0;
}


/**
 * Return a mbr whose first low side is higher than its opposite side -
 * interpreted as an undefined mbr.
 */
 RTREEMBR RTreeNullRect(void)
{
	RTREEMBR mbr;
	int i;

	mbr.bound[0] = (REALTYPE) 1;
	mbr.bound[DIMS_NUMB] = (REALTYPE)-1;
	for (i=1; i<DIMS_NUMB; i++)
		mbr.bound[i] = mbr.bound[i+DIMS_NUMB] = (REALTYPE) 0;
	return mbr;
}

/**
 * Print out the data for a rectangle.
 */
 void RTreePrintRect( RTREEMBR *mbr, int depth)
{
	int i;
	
	_RTreeTabIn(depth);
	fprintf (stdout, "mbr:\n");
	for (i = 0; i < DIMS_NUMB; i++) 
	{
		_RTreeTabIn(depth+1);
		fprintf (stdout, "%f\t%f\n", mbr->bound[i], mbr->bound[i + DIMS_NUMB]);
	}
}

/**
 * Calculate the 2-dimensional area of a rectangle
 */
 REALTYPE RTreeRectArea( RTREEMBR *mbr )
{
	if (INVALID_RECT(mbr))
		return (REALTYPE) 0;

	return (mbr->bound[DIMS_NUMB] - mbr->bound[0]) * (mbr->bound[DIMS_NUMB+1] - mbr->bound[1]);
}


/**
 * Calculate the n-dimensional volume of a rectangle
 */
 REALTYPE RTreeRectVolume( RTREEMBR *mbr )
{
	int i;
	REALTYPE vol = (REALTYPE) 1;

	if (INVALID_RECT(mbr))
		return (REALTYPE) 0;

	for(i=0; i<DIMS_NUMB; i++)
		vol *= (mbr->bound[i+DIMS_NUMB] - mbr->bound[i]);
	assert(vol >= 0.0);
	return vol;
}


/**
 * Precomputed volumes of the unit spheres for the first few dimensions 
 */
const double UnitSphereVolumes[] = {
	0.000000,  /* dimension   0 */
	2.000000,  /* dimension   1 */
	3.141593,  /* dimension   2 */
	4.188790,  /* dimension   3 */
	4.934802,  /* dimension   4 */
	5.263789,  /* dimension   5 */
	5.167713,  /* dimension   6 */
	4.724766,  /* dimension   7 */
	4.058712,  /* dimension   8 */
	3.298509,  /* dimension   9 */
	2.550164,  /* dimension  10 */
	1.884104,  /* dimension  11 */
	1.335263,  /* dimension  12 */
	0.910629,  /* dimension  13 */
	0.599265,  /* dimension  14 */
	0.381443,  /* dimension  15 */
	0.235331,  /* dimension  16 */
	0.140981,  /* dimension  17 */
	0.082146,  /* dimension  18 */
	0.046622,  /* dimension  19 */
	0.025807,  /* dimension  20 */
};

#if DIMS_NUMB > 20
	#error "not enough precomputed sphere volumes"
#endif

#define UnitSphereVolume UnitSphereVolumes[DIMS_NUMB]

/**
 * Calculate the n-dimensional volume of the bounding sphere of a rectangle.
 * The exact volume of the bounding sphere for the given RTREEMBR.
 */
 REALTYPE RTreeRectSphericalVolume( RTREEMBR *mbr )
{
	int i;
	double sum_of_squares=0, radius;

	if (INVALID_RECT(mbr))
		return (REALTYPE) 0;
	
	for (i=0; i<DIMS_NUMB; i++) {
		double half_extent = (mbr->bound[i+DIMS_NUMB] - mbr->bound[i]) / 2;
		sum_of_squares += half_extent * half_extent;
	}
	radius = sqrt(sum_of_squares);
	return (REALTYPE)(pow(radius, DIMS_NUMB) * UnitSphereVolume);
}


/**
 * Calculate the n-dimensional surface area of a rectangle
 */
 REALTYPE RTreeRectSurfaceArea( RTREEMBR *mbr )
{
	int i, j;
	REALTYPE sum = (REALTYPE) 0;

	if (INVALID_RECT(mbr))
		return (REALTYPE) 0;

	for (i=0; i<DIMS_NUMB; i++) {
		REALTYPE face_area = (REALTYPE)1;
		for (j=0; j<DIMS_NUMB; j++)
			/* exclude i extent from product in this dimension */
			if(i != j) {
				REALTYPE j_extent =	mbr->bound[j+DIMS_NUMB] - mbr->bound[j];
				face_area *= j_extent;
			}
			sum += face_area;
	}
	return 2 * sum;
}



/**
 * Combine two rectangles, make one that includes both.
 */
 RTREEMBR RTreeCombineRect( RTREEMBR *rc1, RTREEMBR *rc2 )
{
	int i, j;
	RTREEMBR new_rect;

	assert(rc1 && rc2);

	if (INVALID_RECT(rc1))
		return *rc2;

	if (INVALID_RECT(rc2))
		return *rc1;

	for (i = 0; i < DIMS_NUMB; i++)
	{
		new_rect.bound[i] = MIN(rc1->bound[i], rc2->bound[i]);
		j = i + DIMS_NUMB;
		new_rect.bound[j] = MAX(rc1->bound[j], rc2->bound[j]);
	}
	return new_rect;
}


/**
 * Decide whether two rectangles overlap.
 */
 int RTreeOverlap( RTREEMBR *rc1, RTREEMBR *rc2)
{
	int i, j;
	assert(rc1 && rc2);

	for (i=0; i<DIMS_NUMB; i++)
	{
		j = i + DIMS_NUMB;  /* index for high sides */

		if (rc1->bound[i] > rc2->bound[j] || rc2->bound[i] > rc1->bound[j])
			return FALSE;
	}
	return TRUE;
}


/**
 * Decide whether rectangle r is contained in rectangle s.
 */
 int RTreeContained( RTREEMBR *r, RTREEMBR *s)
{
	int i, j, result;
	assert(r && s);

	/* undefined mbr is contained in any other */
	if (INVALID_RECT(r))
		return TRUE;

	/* no mbr (except an undefined one) is contained in an undef mbr */
	if (INVALID_RECT(s))
		return FALSE;

	result = TRUE;
	for (i = 0; i < DIMS_NUMB; i++)
	{
		j = i + DIMS_NUMB;  /* index for high sides */
		result = result	&& r->bound[i] >= s->bound[i] && r->bound[j] <= s->bound[j];
	}
	return result;
}

/**
 * Split a node.
 * Divides the nodes branches and the extra one between two nodes.
 * Old node is one of the new ones, and one really new one is created.
 * Tries more than one method for choosing a partition, uses best result.
 */
void RTreeSplitNode(HRTREEROOT root, RTREENODE *node, RTREEBRANCH *br, RTREENODE **new_node)
{
	RTREEPARTITION *p;
	int level;

	assert(node && br);
	
	/* load all the branches into a buffer, initialize old node */
	level = node->level;
	_RTreeGetBranches(root, node, br);

	/* find partition */
	p = &root->Partitions[0];

	/* Note: can't use MINFILL(n) below since node was cleared by GetBranches() */
	_RTreeMethodZero(root, p, (level>0 ? MINNODEFILL : MINLEAFFILL));

	/* put branches from buffer into 2 nodes according to chosen partition	*/
	*new_node = RTreeNewNode();
	(*new_node)->level = node->level = level;
	_RTreeLoadNodes(root, node, *new_node, p);

	assert(node->count+(*new_node)->count == p->total);
}


/**
 * Initialize a RTREENODE structure. 
 */
void RTreeInitNode( RTREENODE *node )
{
	int i;
	node->count = 0;
	node->level = -1;
	for (i = 0; i < MAXCARD; i++)
		_RTreeInitBranch(&(node->branch[i]));
}

/** 
 * Make a new node and initialize to have all branch cells empty. 
 */
RTREENODE *RTreeNewNode(void)
{
	RTREENODE *node = (RTREENODE*) malloc(sizeof(RTREENODE));
	assert(node);
	RTreeInitNode(node);
	return node;
}

void RTreeFreeNode( RTREENODE *node )
{
	assert(node);
	free(node);
}


/**
 * Print out the data in a node. 
 */
 void RTreePrintNode( RTREENODE *node, int depth )
{
	int i;
	assert(node);

	_RTreeTabIn(depth);
	fprintf (stdout, "node");

	if (node->level == 0)
		fprintf (stdout, " LEAF");
	else if (node->level > 0)
		fprintf (stdout, " NONLEAF");
	else
		fprintf (stdout, " TYPE=?");

#pragma warning(push)	/* C4311 */
#pragma warning( disable : 4311 )
	fprintf (stdout, "  level=%d  count=%d  address=%o\n", node->level, node->count, (unsigned int) node);

	for (i=0; i<node->count; i++)
	{
		if(node->level == 0) {
			/* _RTreeTabIn(depth); */
			fprintf (stdout, "\t%ld: data = %ld\n", i, (unsigned int)node->branch[i].child);
		}
		else {
			_RTreeTabIn(depth);
			fprintf (stdout, "branch %d\n", i);
			_RTreePrintBranch(&node->branch[i], depth+1);
		}
	}
#pragma warning(pop)
}

/**
 * Find the smallest rectangle that includes all rectangles in branches of a node.
 */
 RTREEMBR RTreeNodeCover( RTREENODE *node )
{
	int i, first_time=1;
	RTREEMBR mbr;
	assert(node);

	RTreeInitRect(&mbr);

	for (i = 0; i < MAXKIDS(node); i++)
	{
		if (node->branch[i].child)
		{
			if (first_time)
			{
				mbr = node->branch[i].mbr;
				first_time = 0;
			}
			else
				mbr = RTreeCombineRect(&mbr, &(node->branch[i].mbr));
		}
	}
	return mbr;
}

/**
 * Pick a branch.  Pick the one that will need the smallest increase
 * in area to accomodate the new rectangle.  This will result in the
 * least total area for the covering rectangles in the current node.
 * In case of a tie, pick the one which was smaller before, to get
 * the best resolution when searching.
 */
 int RTreePickBranch( RTREEMBR *mbr, RTREENODE *node)
{
	RTREEMBR *r;
	int i, first_time = 1;
	REALTYPE increase, bestIncr=(REALTYPE)-1, area, bestArea=0;
	int best=0;
	RTREEMBR tmp_rect;
	assert(mbr && node);

	for (i=0; i<MAXKIDS(node); i++)
	{
		if (node->branch[i].child)
		{
			r = &node->branch[i].mbr;
			area = RTreeRectSphericalVolume(r);
			tmp_rect = RTreeCombineRect(mbr, r);
			increase = RTreeRectSphericalVolume(&tmp_rect) - area;
			if (increase < bestIncr || first_time)
			{
				best = i;
				bestArea = area;
				bestIncr = increase;
				first_time = 0;
			}
			else if (increase == bestIncr && area < bestArea)
			{
				best = i;
				bestArea = area;
				bestIncr = increase;
			}
		}
	}
	return best;
}

/**
 * Add a branch to a node.  Split the node if necessary.
 * Returns 0 if node not split.  Old node updated.
 * Returns 1 if node split, sets *new_node to address of new node.
 * Old node updated, becomes one of two.
 */
 int RTreeAddBranch(HRTREEROOT root, RTREEBRANCH *br, RTREENODE *node, RTREENODE **new_node)
{
	int i;
	assert(br && node);
	
	if (node->count < MAXKIDS(node))  /* split won't be necessary */
	{
		for (i = 0; i < MAXKIDS(node); i++)  /* find empty branch */
		{
			if (node->branch[i].child == NULL)
			{
				node->branch[i] = *br;
				node->count++;
				break;
			}
		}

		return 0;
	}
	
	assert(new_node);
	RTreeSplitNode(root, node, br, new_node);
	
	return 1;
}

/**
 * Disconnect a dependent node. 
 */
 void RTreeDisconnectBranch( RTREENODE *node, int i )
{
	assert(node && i>=0 && i<MAXKIDS(node));
	assert(node->branch[i].child);

	_RTreeInitBranch(&(node->branch[i]));
	node->count--;
}

/**
 * Destroy (free) node recursively. 
 */
 void RTreeDestroyNode ( RTREENODE *node )
{
	int i;

	if (node->level > 0) 
	{
		/* it is not leaf -> destroy childs */
		for ( i = 0; i < NODECARD; i++) 
		{
			if ( node->branch[i].child ) 
				RTreeDestroyNode ( node->branch[i].child );
		}
	}

	/* Free this node */
	RTreeFreeNode( node );
}


/**
 * Create a new rtree index, empty. Consists of a single node. 
 */
HRTREEROOT RTreeCreate(PFNRTreeSearchCallback  pfnSearchCallback)
{
	RTREEROOT * root = (RTREEROOT*) malloc(sizeof(RTREEROOT));
	assert(root);
	root->root_node = RTreeNewNode();
	assert(root->root_node);
	root->root_node->level = 0;		/* leaf */
	root->pfnCallback = pfnSearchCallback;
	return root;
}

/**
 * Destroy a rtree root must be a root of rtree. Free all memory.
 */
void RTreeDestroy(HRTREEROOT root)
{
	RTreeDestroyNode (root->root_node);
	root->root_node = 0;
	free(root);
}

/**
* Search in an index tree for all data rectangles that overlap the argument rectangle.
* Return the number of qualifying data rects.
*/
int RTreeSearch(HRTREEROOT root, RTREEMBR *mbr, void* pfnSearchCallbackParam)
{
	return  _RTreeSearchRect(root->root_node, mbr, root->pfnCallback, pfnSearchCallbackParam);
}


/** 
 * Insert a data rectangle into an index structure.
 * RTreeInsertRect provides for splitting the root;
 * returns 1 if root was split, 0 if it was not.
 * The level argument specifies the number of steps up from the leaf
 * level to insert; e.g. a data rectangle goes in at level = 0.
 * _RTreeInsertRect does the recursion.
 */
int RTreeInsert(HRTREEROOT root, RTREEMBR *mbr, void* tid, int level)
{
#ifdef _DEBUG
	int i;
#endif

	RTREENODE	*newroot;
	RTREENODE	*newnode;
	RTREEBRANCH b;
	
	assert(mbr && root);
	assert(level >= 0 && level <= root->root_node->level);

#ifdef _DEBUG
	for (i=0; i<DIMS_NUMB; i++) 
		assert(mbr->bound[i] <= mbr->bound[DIMS_NUMB+i]);
#endif

	/* root split */
	if (_RTreeInsertRect(root, mbr, tid, root->root_node, &newnode, level))  
	{
		newroot = RTreeNewNode();  /* grow a new root, & tree taller */
		newroot->level = root->root_node->level + 1;
		b.mbr = RTreeNodeCover(root->root_node);
		b.child = root->root_node;
		RTreeAddBranch(root, &b, newroot, NULL);
		b.mbr = RTreeNodeCover(newnode);
		b.child = newnode;
		RTreeAddBranch(root, &b, newroot, NULL);
		root->root_node = newroot;
		
		return 1;
	}

	return 0;
}


/**
 * Delete a data rectangle from an index structure.
 * Pass in a pointer to a RTREEMBR, the tid of the record, ptr to ptr to root node.
 * Returns 1 if record not found, 0 if success.
 * RTreeDeleteRect provides for eliminating the root.
 */
int RTreeDelete(HRTREEROOT root, RTREEMBR *mbr, void* tid)
{
	int		i;
	RTREENODE		*tmp_nptr = NULL;
	RTREELISTNODE	*reInsertList = NULL;
	RTREELISTNODE	*e;

	assert(mbr && root && root->root_node);
	
	if (!_RTreeDeleteRect(mbr, tid, root->root_node, &reInsertList))
	{
		/* found and deleted a data item */

		/* reinsert any branches from eliminated nodes */
		while (reInsertList)
		{
			tmp_nptr = reInsertList->node;

#pragma warning(push)	/* C4311 */
#pragma warning( disable : 4311 )
			for (i = 0; i < MAXKIDS(tmp_nptr); i++)
			{
				if (tmp_nptr->branch[i].child)
				{
					RTreeInsert(root, &(tmp_nptr->branch[i].mbr), (void*)tmp_nptr->branch[i].child, tmp_nptr->level);
				}
			}
#pragma warning(pop)

			e = reInsertList;
			reInsertList = reInsertList->next;
			RTreeFreeNode(e->node);
			_RTreeFreeListNode(e);
		}
		
		/* check for redundant root (not leaf, 1 child) and eliminate */
		if (root->root_node->count == 1 && root->root_node->level > 0)
		{
			for (i = 0; i < NODECARD; i++)
			{
				tmp_nptr = root->root_node->branch[i].child;
				if(tmp_nptr)
					break;
			}
			assert(tmp_nptr);
			RTreeFreeNode(root->root_node);
			root->root_node = tmp_nptr;
		}

		return 0;
	}
	
	return 1;
}