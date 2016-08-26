/****************************************************************************************
 * list.c                                                                               *
 *		Generic sequential linked list node structure -- can hold any type data.        *
 *		cheungmine                                                                      *
 *      Mar. 22, 2008.  All rights reserved.                                            *
 ****************************************************************************************/
#include "list.h"

void 
list_traverse(list_t *in_list, pfunc_list_callback pfcb)
{
	listnode_t* node;
	if (pfcb) {
		node = in_list->head;
		while (node) {
			(*pfcb)(node);
			node = node->next;
		}
	}
}

list_t* 
list_create()
{
	list_t	*list = (list_t*)malloc (sizeof(list_t));
	list->size = 0;
	list->head = list->tail = list->curr = NULL;
	return list;
}

void 
list_destroy(list_t *in_list, pfunc_list_callback  pf)
{
	list_clear(in_list, pf);	
	free(in_list);
}

listnode_t* 
list_node_create(void* data)
{
	listnode_t	*node = (listnode_t*)malloc (sizeof(listnode_t));
	node->next = NULL;
	node->data = data;
	return node;
}

void
list_node_free(listnode_t* node, pfunc_list_callback  pfcb)
{
	listnode_t* next;
	
	while(node){
		next = node->next;
		if (pfcb)
			(*pfcb)(node);
		free(node);
		node = next;
	}
}

listnode_t* 
list_key_create(long key)
{
	listnode_t	*node = (listnode_t*)malloc (sizeof(listnode_t));
	node->next = NULL;
	node->key = key;
	
	return node;
}

listnode_t* 
list_find_prev(const list_t *in_list, const listnode_t *node)
{
	listnode_t* prev;
	assert(node);
	prev = in_list->head;
	if (prev==node)
		return NULL;	
	while(prev && prev->next!=node) {
		prev = prev->next;
	}
	assert(prev);
	return prev;
}

void 
list_push_back(list_t *in_list, listnode_t *node)
{
	node->next = NULL;

	if (in_list->head)
	{
		in_list->tail->next = node;
		in_list->tail = node;
	}
	else
		in_list->head = in_list->tail = node;

	in_list->size++;
}

void 
list_push_front(list_t *in_list, listnode_t *in_node)
{
	in_node->next = in_list->head;
	in_list->head = in_node;

	if (!in_node->next)
		in_list->tail = in_node;

	in_list->size++;
}

void
list_insert_after(list_t *in_list, listnode_t *pos_node, listnode_t *in_node)
{
	assert(pos_node);
	in_node->next = pos_node->next;
	pos_node->next = in_node;

	if (in_list->tail==pos_node)
		in_list->tail = in_node;

	in_list->size++;
}

void
list_insert_before(list_t *in_list, listnode_t *pos_node, listnode_t *in_node)
{
	listnode_t* prev_node;
	prev_node = list_find_prev(in_list, pos_node);
	if (prev_node)
		list_insert_after(in_list, prev_node, in_node);
	else
		list_push_front(in_list, in_node);
}

listnode_t* 
list_pop_front(list_t *in_list)
{
	listnode_t	*pop_node = NULL;
	if (in_list->head)
	{
		pop_node = in_list->head;
		in_list->head = in_list->head->next;
		if (in_list->head == NULL)
			in_list->tail = NULL;
		pop_node->next = NULL;
		in_list->size--;		
	}
	assert(in_list->size >= 0);
	return pop_node;
}

listnode_t* 
list_pop_back(list_t *in_list)
{
	listnode_t	*pop_node = in_list->tail;
	if (in_list->head==pop_node) {
		in_list->size = 0;
		in_list->head = in_list->tail = NULL;
	}
	else {
		assert(pop_node && !pop_node->next);
		in_list->tail = list_find_prev(in_list, pop_node);
		assert(in_list->tail);
		in_list->tail->next = NULL;
		in_list->size--;
	}
	assert(in_list->size >= 0);
	return pop_node;
}

void
list_clear(list_t *in_list, pfunc_list_callback pfcb)
{
	if (in_list->head) {
		listnode_t	*node;
		if (pfcb) {			
			while((node = list_pop_front(in_list))){
				(*pfcb)(node);
				free(node);
			}		
		}
		else {
			while((node = list_pop_front(in_list))){
				free(node);
			}
		}
	}

	assert (!in_list->head && !in_list->tail && in_list->size==0);
}

list_t* 
list_copy(list_t list)
{
	list_t	*newlist = (list_t*)malloc (sizeof(list_t));
	*newlist = list;
	return newlist;
}

void 
list_concat(list_t *first, list_t *second)
{
	if (first->head) {
		if (second->head) {
			first->tail->next = second->head;
			first->tail = second->tail;
			first->size += second->size;
		}
	}
	else
		*first = *second;

	second->head = second->tail = NULL;	
	second->size = 0;
}

size_t 
list_size(const list_t* in_list)
{
	return in_list->size;
}

listnode_t* 
list_node_at(const list_t* in_list, size_t index)
{
	size_t  i=0;
	listnode_t	*node = in_list->head;
	assert(index >=0 && index < in_list->size);
	while (i < index) {
		node = node->next;
		i++;
	}
	return node;
}

listnode_t* 
list_node_erase(list_t* in_list, listnode_t* prev)
{
	listnode_t *erase;
	if (prev==0)
		return list_pop_front(in_list);
	if (prev->next==0)
		return 0;
	assert(prev->next);
	erase = prev->next;
	prev->next = erase->next;
	erase->next = 0;
	in_list->size--;
	return erase;
}

listnode_t*
list_slice(list_t* in_list, listnode_t* begin, listnode_t* end)
{
	listnode_t* prev;
	assert(end);
	prev = list_find_prev(in_list, begin);

	if (prev==NULL) {
		assert(begin==in_list->head);
		in_list->head = end->next;
		end->next = NULL;
		if (in_list->tail == end)
			in_list->tail = in_list->head;
	}
	else {
		prev->next = end->next;
		end->next = NULL;
		if (in_list->tail == end)
			in_list->tail = prev;		
	}
	
	in_list->size = 0;
	prev = in_list->head;
	while(prev) {
		in_list->size++;		
		prev = prev->next;
	}

	return begin;
}

void 
list_reverse(list_t* in_list)
{
	listnode_t* node, *prev, *next;
	
	prev = 0;
	node = in_list->head;
	in_list->tail = node;

	while(node)
	{
		next = node->next;
		node->next = prev;
		prev = node;
		node = next;
	}
	
	in_list->head = prev;
}

/**
 * Rewind list cursor ptr
 */
void 
list_rewind(list_t* in_list) 
{
	in_list->curr = 0;
}

/**
 * Get next node ptr by cursor
 */
listnode_t* list_next_node(list_t* in_list) 
{
	if (in_list->curr == in_list->tail)
		return 0;

	if (in_list->curr == 0)
		in_list->curr = in_list->head;
	else
		in_list->curr = in_list->curr->next;
		
	return in_list->curr;
}
