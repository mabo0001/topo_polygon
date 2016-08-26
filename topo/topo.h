/****************************************************************************** 
 * topo.h - 2D topology implementation
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
 *  Oct. 2000 - First Created  using C++/MFC
 *  Feb. 2008 - Rewritten with C++/STL
 *  Jul. 2008 - Created with C only
 ******************************************************************************/

#ifndef _TOPO_H_INCLUDED
#define _TOPO_H_INCLUDED

#include <stdio.h>
#include <assert.h>

#include "list.h"

/*============================================================================
				               Public Structure
  Abbr:

    T - Type
    V - Vertex
    E - Edge
 
    A - Arc
    P - Polygon

	lst - list_t
	nod - listnode_t

============================================================================*/
#define   TOPO_EPSILON   1.0e-10

/**
 * Hung Arc Flag:
 *   Since it's a unsigned int type, we MUST test with if (flag==-1) or 
 *     if(flag==TOPO_ERRINDEX)
 *     NOT: if(flag<0), it's WRONG! */
#define  TOPO_ERRINDEX  0xffffffff

typedef listnode_t  nodVertexT, nodEdgeT, nodEdge2vT, nodArcT, nodPolygonT;
typedef list_t      lstVertexT, lstEdgeT, lstEdge2vT, lstArcT, lstPolygonT;

/**
 * VertexT* pData = NodData(nodV, VertexT); */
#define NodData(nodePtr,dataType)  ((dataType*)((nodePtr)->data))
#define NodDataC(nodePtr,dataType)  ((const dataType*)((nodePtr)->data))

/**
 * Vertex 2 V structure 
 */
typedef struct                       
{
	double  x, y;
} tVertex2v;

/**
 * Rect 2V structure, must be same as RTREEMBR
 */
typedef struct                       
{
	union{
		struct{
			double  Xmin, Ymin, Xmax, Ymax;
		};
		double mbr[4];
	};  
} tRect2v;

/**
 * Edge2v structure 
 */
typedef struct
{
	tVertex2v  v1;
	tVertex2v  v2;
} tEdge2v;
/**
 * topo vertex structure
 */
typedef struct
{
	union {
		struct{
			double    x;
			double    y;
		};
		tVertex2v    _tV;
	};

	void		*_Id;	   /* vertex id, used by caller, NOT necessary */
	lstEdgeT    *_lstE;	   /* edges share the vertex */
}topo_vertex, VertexT;

/**
 * topo polygon consists of arc list
 * a "C" flaged Polygon consists of counter-clockwise vertices - CCW
 * a "W" flaged Polygon consists of clockwise vertices - CW
 */
typedef struct
{
	void		*_Id;	  /* polygon id, used by caller, NOT necessary */

	lstVertexT  *_lstV;   /* a list of vertex ptr to vertex consists of polygon */

	double       _area;	  /* area of polygon */
	int          _flag;	  /* 0: bad polygon;
						     1: counter-clockwise;
							-1: clockwise */
}topo_polygon, PolygonT;

typedef struct
{
	void		*_Id;	  /* arc id, used by caller, NOT necessary */
	lstVertexT  *_lstE;   /* a list of edges consists of arc */
	double       _length; /* length of arc */
}topo_arc, ArcT;

/**
 * topo edge
 */
typedef struct   
{
	void		*_Id;	  /* edge id, used by caller, NOT necessary */
	
	VertexT     *_V1;	  /* ptr to start vertex */
	VertexT     *_V2;     /* ptr to end vertex */
	nodEdgeT    *_nodEV1; /* ptr to start vertex's edges' list node */
	nodEdgeT    *_nodEV2; /* ptr to end vertex's edges' list node */
	
	PolygonT    *_leftP;  /* ptr to polygon lies left to V1->V2 */
	PolygonT    *_rightP; /* ptr to polygon lies right to V1->V2 */

	ArcT        *_A;	  /* ptr to arc which the edge belongs to */
	nodEdgeT    *_EA;	  /* ptr to arc's edge list node */
}topo_edge, EdgeT;

/**
 * topo net opaque structure
 */
typedef struct _topo_net_t  *topo_net;


/*===========================================================================
							Public Functions
 作者: cheungmine@gmail
 日期: 2008-5
 声明: 拓扑代码, 作者版权所有, 保留所有权利.
===========================================================================*/

void toponet_create(topo_net* tp, double snap /* =0.001? */);
void toponet_free (topo_net t);
void toponet_clear (topo_net t, BOOL clear_edges);
void toponet_set_snap(topo_net t, double snap);
double toponet_get_snap(topo_net t);

void toponet_add_edge (topo_net t, const tVertex2v* v1, const tVertex2v* v2, BOOL split);
void toponet_build (topo_net t, BOOL build_polygons, BOOL build_arcs, BOOL clear_edges);

/* same as:
 * toponet_build(t, TRUE, TRUE, clear_edges); 
 */
void toponet_build_all (topo_net t, BOOL clear_edges);

lstVertexT* toponet_get_vertex_list (topo_net t);
lstEdgeT*  toponet_get_edge_list (topo_net t);
lstArcT*  toponet_get_arc_list (topo_net t);
lstPolygonT*  toponet_get_polygon_list (topo_net t);

void toponet_output_all (topo_net t, FILE* fpOut);
void toponet_output_verts (topo_net t, FILE* fpOut);
void toponet_output_edges (topo_net t, FILE* fpOut);
void toponet_output_arcs (topo_net t, FILE* fpOut);
void toponet_output_polygons (topo_net t, FILE* fpOut);


/*============================================================================*/
#endif  /*_TOPO_H_INCLUDED*/

