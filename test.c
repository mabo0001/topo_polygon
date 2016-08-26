//
// test.c - test topology
// cheungmine@gmail.com
// 2008-3
//
#include <stdio.h>
#include <stdlib.h>

#include "list.h"

void list_use_key()
{
	long k;
	list_t*  lst;
	listnode_t* nod;
	lst = list_create();
	assert(lst);

	for(k=1; k<=10; k++){
		list_push_back(lst, list_key_create(k));
	}

	list_size(lst);
	list_pop_front(lst);
	list_size(lst);
	list_pop_back(lst);
	list_size(lst);

	list_traverse(lst, my_listnode_key_traverse);

	nod = list_slice(lst, lst->head->next->next, list_find_prev(lst, lst->tail));
	list_node_free(nod, 0);

	list_traverse(lst, my_listnode_key_traverse);

	list_destroy(lst, 0);
}

typedef struct 
{
	char * _buff;
	size_t _size;
}blob;

blob* blob_create(size_t n) 
{
	blob* p = (blob*) malloc(sizeof(blob));
	p->_buff = (char*) malloc(n);
	p->_size = n;
	DBG_TRACE("blob_create\n");
	return p;
}

void blob_free(blob* p) 
{
	free(p->_buff);
	p->_size = 0;
	DBG_TRACE("blob_free\n");
}

void my_listnode_blob_free(listnode_t *node)
{
	blob_free((blob*)node->data);
}

void my_listnode_blob_traverse(listnode_t *node)
{
	printf("   _size: %d, _buff: %s\n", ((blob*)node->data)->_size, ((blob*)node->data)->_buff);
}

void list_use_data()
{
	long k;
	list_t*  lst;
	blob* pb;
	listnode_t* nod;
	lst = list_create();
	assert(lst);

	for(k=1; k<=10; k++){
		pb = blob_create(50);
		sprintf_s(pb->_buff, 50, "    this is a blob: %d", k);
		list_push_back(lst, list_node_create(pb));
	}
	list_size(lst);
	
	nod = list_pop_front(lst);
	DBG_TRACE("Use blob here...\n");
	list_node_free(nod, my_listnode_blob_free);
	list_size(lst);

	nod = list_pop_back(lst);
	list_size(lst);
	DBG_TRACE("Use blob here...\n");
	list_node_free(nod, my_listnode_blob_free);
	
	list_traverse(lst, my_listnode_blob_traverse);

	nod = list_slice(lst, lst->head->next->next, list_find_prev(lst, lst->tail));
	list_node_free(nod, my_listnode_blob_free);

	list_traverse(lst, my_listnode_blob_traverse);

	list_destroy(lst, my_listnode_blob_free);
}

int main()
{
	printf("==================list_use_key==================\n");
	list_use_key();

	printf("==================list_use_data==================\n");
	list_use_data();

	printf("-----------------\nPress <Enter> for quit! cheungmine@gmail.com\n");
	getchar();
	return 0;
}