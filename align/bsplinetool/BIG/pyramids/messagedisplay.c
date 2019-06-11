#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"messagedisplay.h"

extern int		MessageDisplay(const char *Message)

{ 
	printf(">>> %s", Message);
	return(0);
}
