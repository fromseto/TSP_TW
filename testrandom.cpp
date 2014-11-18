#include <stdio.h> //rmv coding for randam number access using c++

#include <stdlib.h>

#include <time.h>
#include <iostream>


int main()

{

time_t t;

int i;

srand(time(&t));

for(i=1;i<=10;i++)

std::cout<<(unsigned)rand()%100-90<<"\t";

}