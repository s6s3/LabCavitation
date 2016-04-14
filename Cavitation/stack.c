#include <stdio.h>
#include "stack.h"


int
STACK_push(stack *s, int data) {
	if (s->num >= STACK_SIZE)return -1;
	s->data[s->num++] = data;
	return 0;
}


int 
STACK_pop(stack *s,int *data) {
	if (s->num < 1)return -1;
	*data = s->data[--s->num];
	return 0;
}

int
STACK_length(stack *s) {
	return s->num;
}

void
STACK_init(stack *s) {
	s->num = 0;
}