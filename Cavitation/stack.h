#pragma once

#define STACK_SIZE 100000

typedef struct {
	int num;
	int data[STACK_SIZE];

}stack;

int
STACK_push(stack *s, int data);


int
STACK_pop(stack *s, int *data);

int
STACK_length(stack *s);

void
STACK_init(stack *s);