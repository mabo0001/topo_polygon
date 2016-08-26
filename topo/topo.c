/****************************************************************************** 
 * topo.c - 2D topology implementation
 *
 * Author: cheungmine. cheungmine@gmail.com
 * --------------------------
 *
 * Copyright (C) cheungmine
 *
 * This software is available under the "FreeBSD" license,
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 * --------------------------
 *  Oct. 2000 - First Created  using MFC
 *  Feb. 2008 - Rewritten with C++
 *  Jul. 2008 - Updated with C
 ******************************************************************************/
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <memory.h>

#include "topo.h"

#include "rtree_2d.h"

/*-----------------------------------------------
			  Topo Marcos
-----------------------------------------------*/
#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#ifndef M_HFPI
#define M_HFPI 	1.57079632679489661923
#endif
	
#ifndef M_DBPI
#define M_DBPI 	6.28318530717958647692
#endif

#ifndef MinVal
#define MinVal(a, b)	((a)<(b)? (a):(b))
#endif

#ifndef MaxVal
#define MaxVal(a, b)	((a)>(b)? (a):(b))
#endif

#define  _CreateNod  list_node_create

#define  tSnap		(t->_snap)

/**
 * CCW3V
 *   returns: TRUE is CounterClockWise of (Vi->Vj->Vk); FALSE is for ClockWise
 *     y   Vk
 *      |  /\
 *      | /  \
 *      |/    \    
 *    --Vi====>Vj------> x
 *      |
 */
#define SS3V(Vi, Vj, Vk)    ((Vi)->x*((Vj)->y-(Vk)->y)+(Vj)->x*((Vk)->y-(Vi)->y)+(Vk)->x*((Vi)->y-(Vj)->y))

#define CCW3V(Vi, Vj, Vk)  (SS3V(Vi, Vj, Vk) > 0)

#define  DistSq(v1, v2)   (((v1)->x-(v2)->x)*((v1)->x-(v2)->x) + ((v1)->y-(v2)->y)*((v1)->y-(v2)->y))

#define  DimsV(v)   ((v)->_lstE->size)

#define  EdgeLen(E)   sqrt(DistSq(&E->v1, &E->v2))

/*===========================================================================
							Private Structure/Functions
===========================================================================*/
typedef   pfunc_list_callback   pfnLstNodeDataFree;

/**
 * topo net
 */
typedef struct  _topo_net_t
{
    double		 _snap;	        /* topo tolerance which if ditance between vertices is less than 
						           the vertices will be looked as same one */
	lstVertexT  *_lstV;	 	    /* vertex list */
	lstEdgeT	*_lstE;		    /* edge list */
	lstEdge2vT  *_lstE2v;	    /* edge2v list */
	HRTREEROOT   _hTreeE2v;		/* edge2vs' tree, must create */
	lstEdge2vT  *_lstTmpE2v;	/* temp edge2vs' list */

	lstArcT	    *_lstA;		    /* arc list */
	lstPolygonT *_lstP;		    /* polygon list */

	HRTREEROOT   _hTreeE;		/* edges' tree */
	lstEdgeT    *_lstTmpE;		/* temp edges' list */
}topo_net_t, TopoNetT;

/* get union of two rect */
static void UnionRect(const tRect2v* rc1, const tRect2v* rc2, tRect2v* rcUnion) {
	rcUnion->Xmin = rc1->Xmin<rc2->Xmin? rc1->Xmin : rc2->Xmin;
	rcUnion->Xmax = rc1->Xmax>rc2->Xmax? rc1->Xmax : rc2->Xmax;
	rcUnion->Ymin = rc1->Ymin<rc2->Ymin? rc1->Ymin : rc2->Ymin;  
	rcUnion->Ymax = rc1->Ymax>rc2->Ymax? rc1->Ymax : rc2->Ymax;
}

static int SnapDig(double snap) {
	int n = 0;
	while(snap < 1) {
		n++;
		snap *= 10;
	}
	return n;
}

/*-----------------------------------------------
			 Vertex Private Functions
-----------------------------------------------*/
static VertexT*  VertexCreate(const tVertex2v* tV) {
	VertexT* V = (VertexT*) malloc(sizeof(VertexT));
	assert(V);	
	V->_tV = (*tV);
	V->_lstE = list_create();
	V->_Id = 0;
	return V;
}

static void  VertexFree(VertexT* V) {
	assert(V);
	list_destroy(V->_lstE, 0);
	free(V);
}

/**
 * callback function to free list vertex node data */
static void pfcb_lstV_free(nodVertexT *nodV) {
	VertexFree(NodData(nodV, VertexT));
}

/**
 * get rounding box of vertex */
static void VertexRect(const VertexT* V, double snap, tRect2v* rc) {
	rc->Xmin = V->x - snap;
	rc->Xmax = V->x + snap;
	rc->Ymin = V->y - snap;
	rc->Ymax = V->y + snap;
}

static void  VertexListAdd(topo_net t, nodVertexT *nodV) {
	list_push_back(t->_lstV, nodV);
    NodData(nodV, VertexT)->_Id = (void*)(t->_lstV->size);
}

static BOOL VertexEqual(double snap, const tVertex2v* v1, const tVertex2v* v2) {
	assert(snap >= TOPO_EPSILON);
	if (fabs(v1->x-v2->x) > snap || fabs(v1->y-v2->y) > snap || DistSq(v1, v2) > snap*snap)
		return FALSE;
	return TRUE;
}

/**
 * bind an edge to vertex. Vo is node vertex */
static void VertexBindEdge(VertexT *V, EdgeT *E) {
	EdgeT     *Ei;
	nodEdgeT  *nodEi, *nodEp, *nodEv;
	VertexT   *Vi, *Ve = 0;

	/* Vo->V must be edge E */
	assert(V==E->_V1 || V==E->_V2);
	if (V==E->_V1) {
		nodEv = E->_nodEV1 = _CreateNod(E);
		Ve = E->_V2;
	}
	else {
		nodEv = E->_nodEV2 = _CreateNod(E);
		Ve = E->_V1;
	}
	assert(Ve && Ve != V);

	nodEi = V->_lstE->head;
	nodEp = 0;
	while (nodEi) {
		Ei = NodData(nodEi, EdgeT);
		assert(V==Ei->_V1 || V==Ei->_V2);
		Vi = (V==Ei->_V1? Ei->_V2 : Ei->_V1);
		assert(Vi != V);

		if (CCW3V(V, Vi, Ve)) 
			break;

		nodEp = nodEi;
		nodEi = nodEi->next;
	}
  
	/* create and insert nodE after nodEp */
	if (nodEp)
		list_insert_after(V->_lstE, nodEp, nodEv);
	else
		list_push_front(V->_lstE, nodEv);
}

/*-----------------------------------------------
			  Edge Private Functions
-----------------------------------------------*/
static EdgeT* EdgeCreate(VertexT *V1, VertexT *V2) {
	EdgeT* E = (EdgeT*) malloc(sizeof(EdgeT));
	assert(E);
	
	E->_leftP = E->_rightP = 0;
	E->_V1 = V1;
	E->_V2 = V2;
	E->_nodEV1 = 0;
	E->_nodEV2 = 0;
	E->_Id = 0;
	E->_A = 0;
	E->_EA = 0;
	
	VertexBindEdge(V1, E);
	VertexBindEdge(V2, E);

	return E;
}

static void  EdgeFree(EdgeT* E) {
	assert(E);
	free(E);
}

/**
 * callback function to free list data
 */
static void pfcb_lstE_free(nodEdgeT *nodE) {
	EdgeFree(NodData(nodE, EdgeT));
}

static RTREEMBR* EdgeRect(const EdgeT* E, double snap, tRect2v* rc) {
	rc->Xmin = MinVal(E->_V1->x, E->_V2->x) - snap;
	rc->Xmax = MaxVal(E->_V1->x, E->_V2->x) + snap;
	rc->Ymin = MinVal(E->_V1->y, E->_V2->y) - snap;
	rc->Ymax = MaxVal(E->_V1->y, E->_V2->y) + snap;
	return (RTREEMBR*) rc;
}

static RTREEMBR* Edge2vRect(const tEdge2v* E, double snap, tRect2v* rc) {
	rc->Xmin = MinVal(E->v1.x, E->v2.x) - snap;
	rc->Xmax = MaxVal(E->v1.x, E->v2.x) + snap;
	rc->Ymin = MinVal(E->v1.y, E->v2.y) - snap;
	rc->Ymax = MaxVal(E->v1.y, E->v2.y) + snap;
	return (RTREEMBR*) rc;
}

static void EdgeListAdd(topo_net t, nodEdgeT *nodE) {
	EdgeT* E = NodData(nodE, EdgeT);
	list_push_back(t->_lstE, nodE);
	E->_Id = (void*)t->_lstE->size;

	if (t->_hTreeE) {
		tRect2v   rc;
		RTreeInsert(t->_hTreeE, EdgeRect(E, tSnap, &rc), E, 0);
	}  
}

static int pfcbEdgesSearchCallback(void* data_id, void* pfnParam){
	/** 
	 * data_id is EdgeT*
	 * pfnParam is ptr to t->_lstTmpE */
	list_push_back(pfnParam, _CreateNod(data_id));
	return 1;
}

#define ListCycNextNode(lst, nod)  ((nod)==(lst)->tail? (lst)->head:(nod)->next)
	
static EdgeT *NextEdge(VertexT *Vo, EdgeT *Eo) {
	/**
	 * get edge node joint with Vo */
	nodEdgeT *nodE = (Vo==Eo->_V1?Eo->_nodEV1:Eo->_nodEV2);
	assert(nodE && NodData(nodE, EdgeT)==Eo);
	nodE = ListCycNextNode(Vo->_lstE, nodE);
	assert(nodE);
	return NodData(nodE, EdgeT);
}

static tEdge2v * Edge2vCreate (const tVertex2v* v1, const tVertex2v* v2) {
	tEdge2v* E2v = (tEdge2v*) malloc(sizeof(tEdge2v));
	assert(E2v);
	E2v->v1 = *v1;
	E2v->v2 = *v2;
	return E2v;
}

static void Edge2vFree(tEdge2v* E2v) {
	assert(E2v);
	free (E2v);
}

static void Edge2vListAdd(topo_net t, nodEdge2vT *nodE2v, RTREEMBR* mbr) {
	list_push_back(t->_lstE2v, nodE2v);	
	RTreeInsert(t->_hTreeE2v, mbr, NodData(nodE2v, tEdge2v), 0);
}

static void pfcb_lstE2v_free(nodEdge2vT *nodE2v) {
	Edge2vFree(NodData(nodE2v, tEdge2v));
}

static int SignOf (double v) {
	return (v>0? 1:(v<0?-1:0));
}

/**
 * relation between 2 edges
 * 2 - line a is parallel to line b
 * 1 - edge a has no joint with line b: the start and end vertices of edge a are on the same side of line b
 * 0 - edge a's start or end vertex is on the line b
 *-1 - edge a has an inner joint with line b: the joint vertex is on the edge a except for start and end vertices.
 *     to justify edge a and edge b is joint, must call this function twice as below:
 *              if (Edge2vRelation(a, b)==EDGE_INNERJOINT && Edge2vRelation(b, a)==EDGE_INNERJOINT)
 *                  printf("edge a and edge b is joint!\n");
 */
#define EDGE_PARALLEL		 2	
#define EDGE_NOTJOINT		 1
#define EDGE_VERTEXJOINT	 0
#define EDGE_INNERJOINT		-1

static int Edge2vRelation(const tEdge2v * a, const tEdge2v * b) {
	double	dxb = b->v2.x - b->v1.x;
	double	dyb = b->v2.y - b->v1.y; 
	double	dxa = a->v2.x - a->v1.x;
	double	dya = a->v2.y - a->v1.y;
	
	double	dx1 = a->v1.x - b->v1.x;
	double	dy1 = a->v1.y - b->v1.y;
	double	dx2 = a->v2.x - b->v2.x;	 
	double	dy2 = a->v2.y - b->v2.y;
	
	if ( fabs(dyb*dxa-dxb*dya) < TOPO_EPSILON )
		return  EDGE_PARALLEL;
		
	return (SignOf(dxb*dy1-dyb*dx1) * SignOf(dxb*dy2-dyb*dx2));
}

static BOOL VertexInEdge(const tVertex2v* V, const tEdge2v * E) {
	return ((E->v2.y-V->y)*(V->y-E->v1.y) > 0 || (E->v2.x-V->x)*(V->x-E->v1.x) > 0);
}

/* distance square from vertex to line */
static double DistV2ESq (const tVertex2v* V, const tEdge2v* E, double sqE) {
	double s = SS3V(&E->v1, &E->v2, V);
	return s*s / sqE;
}

/* middle vertex of two vertices
 *   V1======Vmid======V2
 */
#define   MidVertex(V1,V2,Vmid)        (Vmid).x=((V1).x+(V2).x)/2; (Vmid).y=((V1).y+(V2).y)/2


/* get vertex(Vp) inter edge(Ev1--->Ev2) to where dist is sEv1 from Ev1:
 *              
 *  Ev1----sEv1------Vp------>Ev2
 */
static void Edge2vInsVp (const tVertex2v* Ev1, const tVertex2v* Ev2, double sE, double sEv1, tVertex2v* Vp) {
	sEv1 /= sE;
	Vp->x = Ev1->x + (Ev2->x-Ev1->x)*sEv1;
	Vp->y = Ev1->y + (Ev2->y-Ev1->y)*sEv1;
}

/* clips T with E and returns T1 or T2 */
static void Edge2vClip(double sqSnap, tEdge2v* E, const tEdge2v* T, double sqE, double sqT, double sqEv1_T, double sqEv2_T, tEdge2v **pT1, tEdge2v **pT2) {
	double sqEv1_Tv1, sqEv1_Tv2, sqEv2_Tv1, sqEv2_Tv2, sqTv1_F1, sqTv2_F1, sqTv1_F2, sqTv2_F2;
	tVertex2v  F1, F2;
	BOOL  bF1 = FALSE;
	BOOL  bF2 = FALSE;

	if (sqEv1_T <= sqSnap){
		sqEv1_Tv1 = DistSq(&E->v1, &T->v1);
		sqEv1_Tv2 = DistSq(&E->v1, &T->v2);
		sqTv1_F1 = sqEv1_Tv1 - sqEv1_T;
		sqTv2_F1 = sqEv1_Tv2 - sqEv1_T;
		if ((sqTv1_F1+sqTv2_F1)<sqT) {
			double sT = sqrt(sqT);
			double Tv1_F1 = sqrt(sqTv1_F1);
			double Tv2_F1 = sqrt(sqTv2_F1);

			if ((Tv1_F1+Tv2_F1) < (sT+TOPO_EPSILON)) {
				bF1 = TRUE;
				if (Tv1_F1>Tv2_F1)
					Edge2vInsVp(&T->v1, &T->v2, sT, Tv1_F1, &F1);
				else
					Edge2vInsVp(&T->v2, &T->v1, sT, Tv2_F1, &F1);				
			}
		}
	}
	if (sqEv2_T <=sqSnap){
		sqEv2_Tv1 = DistSq(&E->v2, &T->v1);
		sqEv2_Tv2 = DistSq(&E->v2, &T->v2);
		sqTv1_F2 = sqEv2_Tv1 - sqEv2_T;
		sqTv2_F2 = sqEv2_Tv2 - sqEv2_T;
		if ((sqTv1_F2+sqTv2_F2)<sqT) {
			double sT = sqrt(sqT);
			double Tv1_F2 = sqrt(sqTv1_F2);
			double Tv2_F2 = sqrt(sqTv2_F2);

			if ((Tv1_F2+Tv2_F2) < (sT+TOPO_EPSILON)) {
				bF2 = TRUE;
				if (Tv1_F2 > Tv2_F2)
					Edge2vInsVp(&T->v1, &T->v2, sT, Tv1_F2, &F2);
				else
					Edge2vInsVp(&T->v2, &T->v1, sT, Tv2_F2, &F2);				
			}
		}
	}

	if (!bF1 && !bF2)
		return;

	if (bF1 && !bF2) {
		MidVertex(F1, E->v1, E->v1);
		*pT1 = Edge2vCreate(&T->v1, &E->v1);
		*pT2 = Edge2vCreate(&E->v1, &T->v2);		
		return;
	}

	if (bF2 && !bF1) {
		MidVertex(F2, E->v2, E->v2);
		*pT1 = Edge2vCreate(&T->v1, &E->v2);
		*pT2 = Edge2vCreate(&E->v2, &T->v2);
		return;
	}

	assert(bF1 && bF2);
	E->v1 = F1; E->v2 = F2;

	if ((T->v2.y-T->v1.y)*(F2.y-F1.y)>0 || (T->v2.x-T->v1.x)*(F2.x-F1.x)>0) {
		/* Tv1=====F1(Ev1)---E--->(Ev2)F2=====>Tv2 */
		*pT1 = Edge2vCreate(&T->v1, &E->v1);
		*pT2 = Edge2vCreate(&E->v2, &T->v2);
	}
	else {
		/* Tv1=====F2(Ev2)<---E---(Ev1)F1=====>Tv2 */
		*pT1 = Edge2vCreate(&T->v1, &E->v2);
		*pT2 = Edge2vCreate(&E->v1, &T->v2);
	}
}

static void Edge2vSplit(double snap, tEdge2v *T, tEdge2v *E, tEdge2v **pT1, tEdge2v **pT2, tEdge2v **pE1, tEdge2v **pE2) {
	double   sqT, sqE, sqTv1_E, sqTv2_E, sqEv1_T, sqEv2_T;
	double   sqSnap;

	*pT1 = *pT2 = *pE1 = *pE2 = 0;	/* MUST set to ZERO */

	/* check repeat edges */
	if ((VertexEqual(snap,&T->v1,&E->v1)&&VertexEqual(snap,&T->v2,&E->v2)) || (VertexEqual(snap,&T->v1,&E->v2)&&VertexEqual(snap,&T->v2,&E->v1)))
		return;

	/* dist squre */
	sqSnap = snap*snap;
	sqT = DistSq(&T->v1, &T->v2);
	sqE = DistSq(&E->v1, &E->v2);
	
	sqTv1_E = DistV2ESq(&T->v1, E, sqE);
	sqTv2_E = DistV2ESq(&T->v2, E, sqE);
	sqEv1_T = DistV2ESq(&E->v1, T, sqT);
	sqEv2_T = DistV2ESq(&E->v2, T, sqT);

	/* two edges disparts */
	if (MinVal(sqTv1_E, sqTv2_E) > sqSnap && MinVal(sqEv1_T, sqEv2_T) > sqSnap) {
		/* crossing edges */
		if (Edge2vRelation(T, E)==EDGE_INNERJOINT && Edge2vRelation(E, T)==EDGE_INNERJOINT) {
			tVertex2v  Vp;
			tEdge2v    T1, T2, E1, E2;
			double	dx1 = T->v2.x - T->v1.x;
			double	dy1 = T->v2.y - T->v1.y;
			double	dx2 = E->v2.x - E->v1.x;
			double	dy2 = E->v2.y - E->v1.y;		
			double  B = dy1*dx2 - dx1*dy2;		
			assert(fabs(B) >= TOPO_EPSILON); 
			Vp.x = ((E->v1.y-T->v1.y)*dx1*dx2  + T->v1.x*dx2*dy1-E->v1.x*dx1*dy2) / B;
			Vp.y = ((T->v1.x-E->v1.x)*dy1*dy2 - T->v1.y*dx1*dy2+E->v1.y*dy1*dx2) / B;
			T1 = T2 = *T;	/* not using T */
			E1 = E2 = *E;	/* not using E */
			T1.v2 = Vp;	T2.v1 = Vp;	E1.v2 = Vp;	E2.v1 = Vp;
			if (!VertexEqual(snap, &T1.v1, &T1.v2)) *pT1 = Edge2vCreate(&T1.v1, &T1.v2);
			if (!VertexEqual(snap, &T2.v1, &T2.v2)) *pT2 = Edge2vCreate(&T2.v1, &T2.v2);
			if (!VertexEqual(snap, &E1.v1, &E1.v2)) *pE1 = Edge2vCreate(&E1.v1, &E1.v2);
			if (!VertexEqual(snap, &E2.v1, &E2.v2)) *pE2 = Edge2vCreate(&E2.v1, &E2.v2);			
		}
		return;
	}

	/* two verts coincided */
	if (VertexEqual(snap, &T->v1, &E->v1)) {
		if (sqT < sqE) {
			if (sqTv2_E<=sqSnap && VertexInEdge(&T->v2, E)) 
				*pE1 = Edge2vCreate(&T->v2, &E->v2);
		}
		else {
			if (sqEv2_T<=sqSnap && VertexInEdge(&E->v2 ,T)) 
				*pT1 = Edge2vCreate(&E->v2, &T->v2);
		}
		return;
	}	
	if (VertexEqual(snap, &T->v2, &E->v1)) {
		if (sqT < sqE) {
			if (sqTv1_E<=sqSnap && VertexInEdge(&T->v1, E)) 
				*pE1 = Edge2vCreate(&T->v1, &E->v2);
		}
		else {
			if (sqEv2_T<=sqSnap && VertexInEdge(&E->v2 ,T)) 
				*pT1 = Edge2vCreate(&E->v2, &T->v1);
		}
		return;
	}	
	if (VertexEqual(snap, &T->v1, &E->v2)) {
		if (sqT < sqE) {
			if (sqTv2_E<=sqSnap && VertexInEdge(&T->v2, E)) 
				*pE1 = Edge2vCreate(&T->v2, &E->v1);
		}
		else {
			if (sqEv1_T<=sqSnap && VertexInEdge(&E->v1 ,T)) 
				*pT1 = Edge2vCreate(&E->v1, &T->v2);
		}
		return;
	}	
	if (VertexEqual(snap, &T->v2, &E->v2)) {
		if (sqT < sqE) {
			if (sqTv1_E<=sqSnap && VertexInEdge(&T->v1, E))	
				*pE1 = Edge2vCreate(&T->v1, &E->v1);
		}
		else {
			if (sqEv1_T<=sqSnap && VertexInEdge(&E->v1 ,T))	
				*pT1 = Edge2vCreate(&E->v1, &T->v1);
		}
		return;
	}

	/* no verts coincided but closer */
	Edge2vClip(sqSnap, E, T, sqE, sqT, sqEv1_T, sqEv2_T, pT1, pT2);
	Edge2vClip(sqSnap, T, E, sqT, sqE, sqTv1_E, sqTv2_E, pE1, pE2);
}

static BOOL BuildEdge(topo_net t, const tVertex2v* v1, const tVertex2v* v2)
{ 
	EdgeT       *E;
	nodEdgeT	*nodE;
	tRect2v		rc;

	VertexT		*V1 = 0;
	VertexT		*V2 = 0;

	/**
	 * too close of two vertices to add */
	if (VertexEqual(tSnap, v1, v2))
		return FALSE;

	/**
	 * look up the nearest edges */
	if (t->_hTreeE) { /* Have Edges RTree */
		list_clear(t->_lstTmpE, 0);

		VertexRect((const VertexT*)v1, tSnap, &rc);
		RTreeSearch(t->_hTreeE, (RTREEMBR*)&rc, t->_lstTmpE);
		VertexRect((const VertexT*)v2, tSnap, &rc);
		RTreeSearch(t->_hTreeE, (RTREEMBR*)&rc, t->_lstTmpE);

		nodE = t->_lstTmpE->head;
	}
	else
		nodE = t->_lstE->head;

	while(nodE) {
		E = NodData(nodE, EdgeT);

		if (!V1 && VertexEqual(tSnap, (const tVertex2v*) E->_V1, v1))
			V1 = E->_V1;
		if (!V1 && VertexEqual(tSnap, (const tVertex2v*) E->_V2, v1))
			V1 = E->_V2;
		if (!V2 && VertexEqual(tSnap, (const tVertex2v*) E->_V1, v2))
			V2 = E->_V1;
		if (!V2 && VertexEqual(tSnap, (const tVertex2v*) E->_V2, v2))
			V2 = E->_V2;

		if ((V1==E->_V1 && V2==E->_V2) || (V2==E->_V1 && V1==E->_V2))
			return FALSE;
		
		nodE = nodE->next;
	}

	/**
	 * cannot find same one, add new vertex into verts */
	if (!V1) {
		V1 = VertexCreate(v1);
		VertexListAdd(t, _CreateNod(V1));
	}
	if (!V2) {
		V2 = VertexCreate(v2);
		VertexListAdd(t, _CreateNod(V2));
	}
	assert(V1 && V2);

	/* add a new edge */
	E = EdgeCreate(V1, V2);
	EdgeListAdd(t, _CreateNod(E));

	return TRUE;
}

/*-----------------------------------------------
			  Arc Private Functions
-----------------------------------------------*/
static  ArcT* ArcCreate() {
	ArcT* A = (ArcT*)malloc(sizeof(ArcT));
	assert(A);
	A->_Id = 0;
	A->_length = 0;
	A->_lstE = list_create();
	assert(A->_lstE);
	return A;	
}

static void ArcFree(ArcT* A) {
	assert(A);
	list_destroy(A->_lstE, 0);
	free(A);
}

static void pfcb_lstA_free(nodArcT *nodA) {
	ArcFree(NodData(nodA, ArcT));
}

static void ArcAddEdge(ArcT *A, EdgeT *E) {
	E->_EA = _CreateNod(E);
	list_push_back(A->_lstE, E->_EA);
	E->_A = A;	
}

static double ArcLength(ArcT *A) {
	nodEdgeT   *nodE;
	EdgeT      *E;
	A->_length = 0;
	nodE = A->_lstE->head;
	while (nodE) {
		E = NodData(nodE, EdgeT);
		A->_length += sqrt(DistSq(E->_V1, E->_V2));
		nodE = nodE->next;
	}
	return A->_length;
}

static void ArcListAdd (topo_net t, ArcT *A) {
	list_push_back(t->_lstA, _CreateNod(A));
	A->_Id = (void*) t->_lstA->size;
	ArcLength(A);
}

/*-----------------------------------------------
			  Polygon Private Functions
-----------------------------------------------*/
static PolygonT* PolygonCreate() {
	PolygonT* P = (PolygonT*)malloc(sizeof(PolygonT));
	assert(P);
	P->_area = 0;
	P->_flag = 0;
	P->_lstV = list_create();
	assert(P->_lstV);
	return P;
}

static void PolygonFree(PolygonT* P) {
	assert(P);	
	/**
	 * since we do NOT call VertexCreate for vertex node, 
	 * so we need NOT set a vertex free callback */
	list_destroy(P->_lstV, 0);
	free(P);
}

static void PolygonAddVertex(PolygonT *P, VertexT *V) {
	list_push_back(P->_lstV, _CreateNod(V));
}

static void PolygonReverse(PolygonT *P) {	
	list_reverse(P->_lstV);
	P->_flag *= (-1);	
}

static void pfcb_lstP_free(nodPolygonT *nodP) {
	PolygonFree(NodData(nodP, PolygonT));
}

static void PolygonBuild (PolygonT *P, EdgeT *E0, VertexT *V0)
{
	EdgeT    *E = E0;
	VertexT  *V = V0;
	assert(E && V);

	PolygonAddVertex(P, V);

	do{
		if (E->_V1==V) {
			assert(E->_leftP==0);
			E->_leftP = P;
			V = E->_V2;
		}
		else {
			assert(E->_rightP==0);
			E->_rightP = P;
			V = E->_V1;
		}

		assert(V);
		PolygonAddVertex(P, V);

		E = NextEdge(V, E);
		assert(E);

	}while (((V0==E->_V1&&E->_leftP==0)||(V0==E->_V2&&E->_rightP==0)) || (V != V0));
}

static double PolygonArea(PolygonT* P) {
	nodVertexT   *nodV, *nodPrev, *nodNext;
	VertexT      *V, *Prev, *Next;

	P->_area = 0;
	P->_flag = 0;
	
	if (P->_lstV->size > 3) {
		nodPrev = P->_lstV->head;
		nodV = nodPrev->next;
		nodNext = nodV->next;
		assert(nodPrev && nodV && nodNext);

		while(nodNext) {
			Prev = NodData(nodPrev, VertexT);
			V = NodData(nodV, VertexT);
			Next = NodData(nodNext, VertexT);
			P->_area += V->x*(Next->y - Prev->y);
			nodPrev = nodV;
			nodV = nodNext;
			nodNext = nodNext->next;
		}

		assert(!nodNext);
		nodNext = P->_lstV->head->next;
		assert(nodV==P->_lstV->tail);

		Prev = NodData(nodPrev, VertexT);
		V = NodData(nodV, VertexT);
		Next = NodData(nodNext, VertexT);
		P->_area += V->x*(Next->y - Prev->y);
	
		P->_flag = (P->_area>0? 1 : (P->_area<0? -1 : 0));
		P->_area = fabs(P->_area)/2;
	}
	return P->_area;
}

static void PolygonListAdd(topo_net t, PolygonT *P) {
	list_push_back(t->_lstP, _CreateNod(P));
	P->_Id = (void*)t->_lstP->size;
	PolygonArea(P);
}

/*-----------------------------------------------
			  Toponet Private Functions
-----------------------------------------------*/
static void ToponetBuildEdges (topo_net t)
{
	nodEdge2vT  *nodE2v;
	tEdge2v		*E2v;

	assert(t->_lstE->size == 0);

	nodE2v = t->_lstE2v->head;
	while (nodE2v) {
		E2v = NodData(nodE2v, tEdge2v);
		BuildEdge(t, &E2v->v1, &E2v->v2);
		nodE2v = nodE2v->next;
	}
}

static void ToponetBuildArcs(topo_net t)
{
#define NextV(E,V)   (E->_V1==V?E->_V2:E->_V1)
#define NextE(V,E)   (E==NodData(V->_lstE->head,EdgeT)?NodData(V->_lstE->tail,EdgeT):NodData(V->_lstE->head,EdgeT))
	nodVertexT   *nodV;
	nodEdgeT     *nodE;
	EdgeT        *E;
	VertexT		 *V, *V0;
			
	assert(t->_lstA->size == 0);

	/* build node arcs or hung arcs */
	nodV = t->_lstV->head;
	while (nodV) {
		V0 = NodData(nodV, VertexT);
		
		if (DimsV(V0) != 2) {
			nodE = V0->_lstE->head;

			while (nodE) {
				E = NodData(nodE, EdgeT);
				
				if (E->_A==0){
					/* create new arc */
					ArcT *A = ArcCreate();

					ArcAddEdge(A, E);
					V = NextV(E,V0);
					
					while (DimsV(V)==2) {
						E = NextE(V,E);
						ArcAddEdge(A, E);
						V = NextV(E,V);
					}

					ArcListAdd(t, A);
				}
				
				nodE = nodE->next;
			}			
		}
		nodV = nodV->next;
	}

	/* build self-closed arcs MUST after building node arcs */
	nodV = t->_lstV->head;
	while (nodV) {
		V0 = NodData(nodV, VertexT);

		if (DimsV(V0) == 2) {
			E = NodData(V0->_lstE->head, EdgeT);
			if (E->_A == 0) {
				ArcT *A = ArcCreate();
				ArcAddEdge(A, E);
				V = NextV(E,V0);

				while (DimsV(V)==2 && V != V0) {
					E = NextE(V,E);
					ArcAddEdge(A, E);
					V = NextV(E,V);
				}

				ArcListAdd(t, A);
			}
		}
		nodV = nodV->next;
	}

#undef NextE
#undef NextV
}

static void ToponetBuildPolygons(topo_net t)
{
	PolygonT    *P;
	nodEdgeT	*nodE;
	EdgeT       *E;

	assert(t->_lstP->size==0);

	nodE = t->_lstE->head;
	while (nodE) {
		E = NodData(nodE, EdgeT);

		if (E->_leftP==0) {
			P = PolygonCreate();
			PolygonBuild(P, E, E->_V1);
			assert(E->_leftP);
			PolygonListAdd(t, P);
		}

		assert(E==NodData(nodE, EdgeT));

		if (E->_rightP==0) {
			P = PolygonCreate();
			PolygonBuild(P, E, E->_V2);
			assert(E->_rightP);
			PolygonListAdd(t, P);
		}

		nodE = nodE->next;
	}
}
/*===========================================================================
							Public Algorithm Functions
===========================================================================*/

void toponet_create(topo_net* pTN, double snap)
{
	topo_net t = (topo_net) calloc(1, sizeof(TopoNetT));
	assert(t);

	toponet_set_snap(t, snap);
	
	t->_lstV = list_create(); assert(t->_lstV);
	t->_lstE = list_create(); assert(t->_lstE);
	t->_lstE2v = list_create(); assert(t->_lstE2v);
	t->_lstA = list_create(); assert(t->_lstA);
	t->_lstP = list_create(); assert(t->_lstP);
	t->_lstTmpE = list_create(); assert(t->_lstTmpE);
	t->_lstTmpE2v = list_create(); assert(t->_lstTmpE2v);

	t->_hTreeE = 0;
	
	t->_hTreeE2v = RTreeCreate(pfcbEdgesSearchCallback);
	assert(t->_hTreeE2v);
	
	*pTN = t;
}

void toponet_free(topo_net t)
{
	assert(t);

	list_destroy(t->_lstV, pfcb_lstV_free);
	list_destroy(t->_lstE, pfcb_lstE_free);
	list_destroy(t->_lstE2v, pfcb_lstE2v_free);
	list_destroy(t->_lstA, pfcb_lstA_free);
	list_destroy(t->_lstP, pfcb_lstP_free);	
	list_destroy(t->_lstTmpE, 0);
	list_destroy(t->_lstTmpE2v, 0);

	if (t->_hTreeE) {
		RTreeDestroy(t->_hTreeE);
		t->_hTreeE = 0;
	}
	if (t->_hTreeE2v) {
		RTreeDestroy(t->_hTreeE2v);
		t->_hTreeE2v = 0;
	}
	free(t);
}

void toponet_set_snap(topo_net t, double snap)
{
	t->_snap = (snap < TOPO_EPSILON? TOPO_EPSILON : snap);
}

double toponet_get_snap(topo_net t)
{
	return t->_snap;
}

/**
 * this is the most important method for splitting and cleaning edges
 */
void toponet_add_edge(topo_net t, const tVertex2v* v1, const tVertex2v* v2, BOOL split)
	{
		lstEdge2vT *lstE;
		tEdge2v    *E, *T, *E1, *E2, *T1, *T2;
		tRect2v     rc;

		nodEdge2vT *nodT;
		nodEdge2vT *nodE, *nodEprev;
		BOOL		bBreakT;
	
	/* must first check: EDGE IS TOO SMALL */
	if (VertexEqual(tSnap, v1, v2))
		return;	

	/* create new edge */
	E = Edge2vCreate(v1, v2);
	assert(E);

	if (!split) {
		Edge2vListAdd(t, _CreateNod(E), Edge2vRect(E, tSnap, &rc));
		return;
	}
	
	/* search near edges in edges tree */
	assert(t->_hTreeE2v);
	list_clear(t->_lstTmpE2v, 0);	/* must clear TEMP list first! */
	RTreeSearch(t->_hTreeE2v, Edge2vRect(E, tSnap, &rc), t->_lstTmpE2v);

	if (t->_lstTmpE2v->size==0) {
		Edge2vListAdd(t, _CreateNod(E), (RTREEMBR*)&rc);
		return;
	}

	/* create new edge list lstE whose contains edges to be added in */
	assert(t->_lstTmpE2v->size>0);
	lstE = list_create();
	assert(lstE);	
	list_push_back(lstE, _CreateNod(E));
	
	nodT = t->_lstTmpE2v->head;
	while (nodT) {
		T = NodData(nodT, tEdge2v);
		bBreakT = FALSE;

		nodE = lstE->head;
		nodEprev = 0;
		while(nodE && T) {
			E = NodData(nodE, tEdge2v);
			
			bBreakT = FALSE;
			Edge2vSplit(tSnap, T, E, &T1, &T2, &E1, &E2);

			/* if T is split into T1, T2, add to front of _lstE2v and drop T at once */
			if (T1 || T2) {
				if (T1) Edge2vListAdd(t, _CreateNod(T1), Edge2vRect(T1, tSnap, &rc));
				if (T2) Edge2vListAdd(t, _CreateNod(T2), Edge2vRect(T2, tSnap, &rc));

				RTreeDelete(t->_hTreeE2v, Edge2vRect(T, tSnap, &rc), T);
				T->v1 = T->v2;	/* mark invalid edge by set v1=v2 since we cannot get actual edge node's addr */

				/* move next T? */
				bBreakT = TRUE;
				nodT = nodT->next;
				T = (nodT? NodData(nodT, tEdge2v) : 0);				
			}
									
			/* if E is split into E1, E2, add to back of lstE and drop E */
			if (E1 || E2) {
				if (E1) list_push_back(lstE, _CreateNod(E1));
				if (E2) list_push_back(lstE, _CreateNod(E2));

				nodE = nodE->next;
				list_node_free(list_node_erase(lstE, nodEprev), pfcb_lstE2v_free);
			}
			else {
				assert(!E1 && !E2);
				nodEprev = nodE;
				nodE = nodE->next;
			}
		}
	
		/* traverse to next node */
		if (!nodT)
			break;

		if (!bBreakT)
			nodT = nodT->next;
	}

	/* clear TEMP list after use */
	list_clear(t->_lstTmpE2v, 0);

	/* add new edges into edges list */
	list_concat(t->_lstE2v, lstE);
	list_destroy(lstE, pfcb_lstE2v_free);
}

void toponet_clear (topo_net t, BOOL clear_edges)
{
	list_clear(t->_lstV, pfcb_lstV_free);
	list_clear(t->_lstE, pfcb_lstE_free);
	list_clear(t->_lstA, pfcb_lstA_free);	
	list_clear(t->_lstP, pfcb_lstP_free);

	if (t->_hTreeE) RTreeDestroy(t->_hTreeE);
	t->_hTreeE = RTreeCreate(pfcbEdgesSearchCallback);
	assert(t->_hTreeE);

	if (clear_edges) {
		RTreeDestroy(t->_hTreeE2v);
		t->_hTreeE2v = RTreeCreate(pfcbEdgesSearchCallback);
		list_clear(t->_lstE2v, pfcb_lstE2v_free);
	}
}

void toponet_build (topo_net t, BOOL build_polygons, BOOL build_arcs, BOOL clear_edges)
{
	toponet_clear (t, FALSE);

	ToponetBuildEdges(t);

	if (build_arcs)
		ToponetBuildArcs(t);

	if (build_polygons)
		ToponetBuildPolygons(t);

	list_clear(t->_lstTmpE, 0);
	list_clear(t->_lstTmpE2v, 0);

	if (clear_edges) {
		RTreeDestroy(t->_hTreeE2v);
		t->_hTreeE2v = RTreeCreate(pfcbEdgesSearchCallback);
		list_clear(t->_lstE2v, pfcb_lstE2v_free);
	}
}

void toponet_build_all (topo_net t, BOOL clear_edges)
{
	toponet_build (t, TRUE, TRUE, clear_edges);
}
/*===========================================================================
							Public Output Functions
===========================================================================*/
#define  AddrOf(ptr)   ((size_t)(void*)(ptr))

lstVertexT* toponet_get_vertex_list (topo_net t){ return t->_lstV; }
lstEdgeT*  toponet_get_edge_list (topo_net t){ return t->_lstE; }
lstArcT*  toponet_get_arc_list (topo_net t){ return t->_lstA; }
lstPolygonT*  toponet_get_polygon_list (topo_net t){ return t->_lstP; }

void toponet_output_all (topo_net t, FILE* fp)
{
	fprintf(fp, "********* TOPONET OUTPUT - cheungmine@gmail.com *********\n");
	
	toponet_output_verts (t, fp);
	toponet_output_edges (t, fp);
	toponet_output_arcs (t, fp);
	toponet_output_polygons (t, fp);
}

void toponet_output_verts (topo_net t, FILE* fp)
{
	const nodVertexT   *nodV;
	const VertexT      *V;
	const nodEdgeT     *nodE;
	const EdgeT        *E;
	int   dig = SnapDig(tSnap);
  
	fprintf(fp, "*Vert: %d\n  V         X          Y     E{n}\n", t->_lstV->size);
	fprintf(fp, "----------------------------------------------------------------\n");
  
	nodV = t->_lstV->head;
	while(nodV) {
		V = NodDataC(nodV, VertexT);
		
		fprintf(fp, "  %d  %10.*f %10.*f    {%d}", AddrOf(V->_Id), dig, V->x, dig, V->y, DimsV(V));

		nodE = V->_lstE->head;
		while(nodE) {
			E = NodDataC(nodE, EdgeT);			
			fprintf(fp, "=[%d]", AddrOf(E->_Id));
			nodE = nodE->next;
		}		
		fprintf(fp, "\n");

		nodV = nodV->next;
	}
	fprintf(fp, "----------------------------------------------------------------\n\n");
}

void toponet_output_edges (topo_net t, FILE* fp)
{
	const nodEdgeT     *nodE;
	const EdgeT        *E;
  
	fprintf(fp, "*Edge: %d\n  E    V1    V2    LeftP  RightP   A\n", t->_lstE->size);
	fprintf(fp, "----------------------------------------------------------------\n");
  
	nodE = t->_lstE->head;
	while(nodE) {
		E = NodDataC(nodE, EdgeT);

		fprintf(fp, " [%d]    %d     %d     <%d>    <%d>    (%d)\n", AddrOf(E->_Id), AddrOf(E->_V1->_Id), AddrOf(E->_V2->_Id), AddrOf(E->_leftP->_Id), AddrOf(E->_rightP->_Id), AddrOf(E->_A->_Id));

		nodE = nodE->next;
	}
	fprintf(fp, "----------------------------------------------------------------\n\n");
}

void toponet_output_arcs (topo_net t, FILE* fp)
{
	const nodArcT     *nodA;
	const ArcT        *A;
	const nodEdgeT    *nodE;
	const EdgeT       *E;
	int   dig = SnapDig(tSnap);
  
	fprintf(fp, "*Arc: %d\n  A         Length    E{n}\n", t->_lstA->size);
	fprintf(fp, "----------------------------------------------------------------\n");
  
	nodA = t->_lstA->head;
	while(nodA) {
		A = NodDataC(nodA, ArcT);

		fprintf(fp, " (%d)    %10.*f     {%d}", AddrOf(A->_Id), dig, A->_length, A->_lstE->size);
		nodE = A->_lstE->head;
		while(nodE) {
			E = NodDataC(nodE, EdgeT);			
			fprintf(fp, "=[%d]", AddrOf(E->_Id));
			nodE = nodE->next;
		}		
		fprintf(fp, "\n");

		nodA = nodA->next;
	}
	fprintf(fp, "----------------------------------------------------------------\n\n");
}

void toponet_output_polygons (topo_net t, FILE* fp)
{
	const nodPolygonT  *nodP;
	const PolygonT     *P;
	const nodVertexT   *nodV;
	const VertexT      *V;
	int   dig = SnapDig(tSnap);
  
	fprintf(fp, "*Polygon: %d\n  P        Area    Flag  V{n}\n", t->_lstP->size);
	fprintf(fp, "----------------------------------------------------------------\n");
  
	nodP = t->_lstP->head;
	while(nodP) {
		P = NodDataC(nodP, PolygonT);

		fprintf(fp, " <%d>  %10.*f   %2d     {%d}", AddrOf(P->_Id), dig, P->_area, P->_flag, P->_lstV->size);

		nodV = P->_lstV->head;
		while(nodV) {
			V = NodDataC(nodV, VertexT);			
			fprintf(fp, "-%d", AddrOf(V->_Id));
			nodV = nodV->next;
		}		
		fprintf(fp, "\n");
		
		nodP = nodP->next;
	}
	fprintf(fp, "----------------------------------------------------------------\n\n");
}

