#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Delcare struct room and room name 
struct room {
    int length; 
    int width; 
}bathroom1, bedroom1, bedroom2_1, bedroom2_2, 
kitchen1, kitchen2, living1, living2;

//Return name and value of larger house using struct 
struct comp{
    char name[20]; 
}p[2];


int main(){
    int spr1, spr2;
    //House 1
    printf("Enter House 1: \n");
    printf("Enter bathroom, bedroom, kitchen, living length: ");
    scanf("%d %d %d %d", &bathroom1.length, &bedroom1.length, &kitchen1.length, &living1.length);
    printf("Enter bathroom, bedroom, kitchen, living width: ");
    scanf("%d %d %d %d", &bathroom1.width, &bedroom1.width, &kitchen1.width, &living1.width);
    
    //House 2
    printf("Enter House 2: \n");
    printf("Enter bedroom 1, bedroom 2, kitchen, living length: ");
    scanf("%d %d %d %d", &bedroom2_1.length, &bedroom2_2.length, &kitchen2.length, &living2.length);
    printf("Enter bedroom 1, bedroom 2, kitchen, living width: ");
    scanf("%d %d %d %d", &bedroom2_1.width, &bedroom2_2.width, &kitchen2.width, &living2.width);

    // Calculate Area 1 length and width
    spr1 =(bathroom1.length+ bathroom1.width)+(bedroom1.length*bedroom1.width)
    +(living1.length*living1.width)+(kitchen1.length*kitchen1.width);

    // Calculate Area 2 length and width
    spr2 =(bedroom2_1.length+ bedroom2_1.width)+(bedroom2_2.length*bedroom2_2.width)
    +(living1.length*living1.width)+(kitchen1.length*kitchen1.width);

    printf("The House 1 square feet is: %d \nThe House 2 square feet is: %d", spr1, spr2);
    printf("\n");

    strcpy(p[0].name,"HOUSE 1");
    strcpy(p[1].name,"HOUSE 2");

    if (spr1 > spr2){
        printf ("The larger house is: %s", p[0].name);
    }
    else
    {
        printf("The larger house is: %s", p[1].name);
    }
    return 0;
}   

