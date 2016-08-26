//
// test.c - test topology
// cheungmine@gmail.com
// 2008-3
//

/////////////////////////////////////
// �ڴ�й©���
// ����Ҫ���ĵط��������:
//        _CrtDumpMemoryLeaks();
// ����3��Ĵ����ܸı�
#define _CRTDBG_MAP_ALLOC
#include <stdio.h>

#include<stdlib.h>
#include <share.h>
//#include<crtdbg.h>

#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include "../topo.h"

#include <windows.h>

void  TestSampleA() 
{
#define NUMEDGES 1000

	const tEdge2v  edges[] = {
		{{50,50},{100,0}},
		{{0,0},{100,0}},
		{{0,0},{0,100}},
		{{0,0},{-100,0}},
		{{0,0},{-50,50}},
		{{-50,50},{0,100}},
		{{-50,50},{-100,0}},
		{{200,200},{300,200}},
		{{300,300},{300,200}},
		{{200,200},{300,300}}	
	};

	DWORD  t0;
	int  i, j;
	topo_net  topo;
	BOOL      useRTree = TRUE;

	int N = sizeof(edges)/sizeof(edges[0]);
	
	toponet_create(&topo, 0.001);
		
	printf("TOPO: add %d edges...", NUMEDGES*N);
	t0 = GetTickCount();
	
	// test for big data
	for (j=0; j<NUMEDGES; j++){
		for (i=0; i<N; i++)
		{
			tVertex2v  v1 = edges[i].v1;
			tVertex2v  v2 = edges[i].v2;

			v1.x += j*100;
			v1.y += j*100;
			
			v2.x += j*100;
			v2.y += j*100;

			toponet_add_edge(topo, &v1, &v2, TRUE);
		}
	}

	printf("%d msec.\n\nTOPO: build...", GetTickCount()-t0);
	t0 = GetTickCount();
	toponet_build_all(topo, TRUE);

	printf("%d msec.\n\nTOPO: build OK, press <Enter> for output\n", GetTickCount()-t0);
	getchar();

	toponet_output_verts (topo, stdout);
	toponet_output_edges (topo, stdout);
	toponet_output_arcs (topo, stdout);
	toponet_output_polygons (topo, stdout);

	toponet_free(topo);

	//_CrtDumpMemoryLeaks();
    printf("Press <Enter> for Quit\n");
	getchar();	
}


#define  TEST_FILE  "C:\\topo.lin"

/*void TestSampleB()
{
	FILE  *fp = 0;
	FILE  *fpOut = 0;
	topo_net  topo;
	tVertex2v v1, v2;
	int   i, nEdges = 0;
	float fSnap = 0.0f;

	fopen_s(&fp, TEST_FILE, "r");
	assert(fp);
	
	fscanf_s(fp, "%d,%f", &nEdges, &fSnap);

	toponet_create(&topo, fSnap);
		
	printf("TOPO: add %d edges...\n", nEdges);
		
	for (i=0; i<nEdges; i++){
		fscanf_s(fp, "%lf,%lf,%lf,%lf\n", &v1.x, &v1.y, &v2.x, &v2.y);
		toponet_add_edge (topo, &v1, &v2, 1);
	}
	fclose(fp);

	toponet_build_all(topo, TRUE);

	fopen_s(&fpOut, "C:\\topo_output.txt", "w");
	assert(fpOut);

	toponet_output_all (topo, fpOut);
	fclose(fpOut);

	toponet_free(topo);

	//_CrtDumpMemoryLeaks();
    printf("Press <Enter> for Quit\n");
	getchar();	
}*/

void TestSampleC()
{
	/* 
	          y
	       |     |
	 (2)---+-----+----
           |     |
	 (1)---+-----+---- ===>x
	       |     |
          (3)   (4)
	 */
	const tEdge2v  edges2[] = {
		{{-190,0},{190,0}},	    // (0) error edge
		{{-80,0.0008},{100,-0.0002}},		// (1)
		{{-100,100},{100,100}}, // (2)
		{{-50,-50},{-50,200}},  // (3)
		{{50,-50},{50,200}}     // (4)
	};

	const tEdge2v  edges[] = {
		{{100,0},{0,0}},
		{{100,0},{100,100}},
		{{0,0},{100,100}}	
	};

	int       i, j=0;
	topo_net  topo;
	int N = sizeof(edges)/sizeof(edges[0]);
	
	toponet_create(&topo, 0.001);
		
	for (i=0; i<N; i++)
		toponet_add_edge (topo, &edges[i].v1, &edges[i].v2, TRUE);

	toponet_build_all (topo, TRUE);
	toponet_output_all (topo, stdout);
	toponet_free (topo);

	//_CrtDumpMemoryLeaks();

    printf("Press <Enter> for Quit\n");
	getchar();	
}


void TestSampleD()
{
/*
                {1,4} {2,4} {3,4}  
        {0,4} +----+----+----+---+ {4,4}
        	  |         |        |
	    {0,3} +         +{2,3}   + {4,3}
	          |  {1,2}{2,2}{3,2} |
        {0,2} +----+----+----+---+ {4,2}
	          |         |        |
	    {0,1} +         +{2,1}   + {4,1}
	          |         |        |
        {0,0} +----+----+----+---+ {4,0}
                 {1,0} {2,0} {3,0}
	 */
	
	const tEdge2v  edges[] = {
		{ {0,0}, {0,1}},
		{ {0,1}, {0,2}},
		{ {0,2}, {0,3}},
		{ {0,3}, {0,4}},
		{ {0,4}, {1,4}},
		{ {1,4}, {2,4}},
		{ {2,4}, {3,4}},
		{ {3,4}, {4,4}},
		{ {4,4}, {4,3}},
		{ {4,3}, {4,2}},
		{ {4,2}, {4,1}},
		{ {4,1}, {4,0}},
		{ {4,0}, {3,0}},
		{ {3,0}, {2,0}},
		{ {2,0}, {1,0}},
		{ {1,0}, {0,0}},

		{ {2,2}, {2,3}},
		{ {2,3}, {2,4}},
		{ {2,2}, {3,2}},
		{ {3,2}, {4,2}},
		{ {2,2}, {2,1}},
		{ {2,1}, {2,0}},
		{ {2,2}, {1,2}},
		{ {1,2}, {0,2}}
	};

	const tEdge2v  edges2[] = {
		{ {0,0}, {4,0}},
		{ {0,2}, {4,2}},
		{ {0,4}, {4,4}},

		{ {0,0}, {0,4}},
		{ {2,0}, {2,4}},
		{ {4,0}, {4,4}}
	};

	int       i, j=0;
	topo_net  topo;
	int N = sizeof(edges)/sizeof(edges[0]);
	
	toponet_create(&topo, 0.001);
		
	for (i=0; i<N; i++)
		toponet_add_edge (topo, &edges[i].v1, &edges[i].v2, TRUE);

	toponet_build_all (topo, TRUE);
	toponet_output_all (topo, stdout);
	toponet_free (topo);

	//_CrtDumpMemoryLeaks();

    printf("Press <Enter> for Quit\n");
	getchar();	
}

int main()
{
	// TestSampleA();
	// TestSampleB();
	// TestSampleC();

	TestSampleD();

	return 0;
}