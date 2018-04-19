#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ConfigParser.h"

int config_get_int(char *search_section, char *search_key){
	int i;
	for(i=0; i<config_num; i++){
		if (!strcmp(search_section, config_section[i]) && !strcmp(search_key, config_key[i])){
			return atoi(config_params[i]);
		}
	}
	printf("There is no such combination of section (%s) and key(%s).\n",search_section,search_key);
	exit(1);
}

float config_get_float(char *search_section, char *search_key){
	int i;
	for(i=0; i<config_num; i++){
		if (!strcmp(search_section, config_section[i]) && !strcmp(search_key, config_key[i])){
			return atof(config_params[i]);
		}
	}
	printf("There is no such combination of section (%s) and key(%s).\n",search_section,search_key);
	exit(1);
}

char *config_get_string(char *search_section, char *search_key){
	int i;
	for(i=0; i<config_num; i++){
		if (!strcmp(search_section, config_section[i]) && !strcmp(search_key, config_key[i])){
			return config_params[i];
		}
	}
	printf("There is no such combination of section (%s) and key(%s).\n",search_section,search_key);
	exit(1);
}

int *config_get_intarr(char *search_section, char *search_key, int *ssize){
	int i,j;
	int size=0;
	for(i=0; i<config_num; i++){
		if (!strcmp(search_section, config_section[i]) && !strcmp(search_key, config_key[i])){
			char *token = strtok(config_params[i], " ,\n");
			int *arr = malloc(sizeof(int));
			while(token){
				if (arr==NULL){
					printf("Failed while allocating memory for int array\n");
					exit(1);
				}
				arr[size] = atoi(token);
				token = strtok(NULL, " ,\n");
				size++;
				arr = realloc(arr, (size+1)*sizeof(int));
			}
			*ssize = size;
			return arr;
		}
	}
	printf("There is no such combination of section (%s) and key (%s).\n",search_section,search_key);
	exit(1);
}

void read_config(char *filename){
	FILE *fp = fopen(filename,"r");
	if (fp==NULL){
		printf("%s cannot be opened.\n",filename);
		exit(1);
	}
	char * line = NULL;
	char * temp_section = malloc(1024 * sizeof(char));;
	char * temp_key = NULL;
	char * temp_params = NULL;
	ssize_t lread;
	size_t len = 0;
	config_num=0;
	
	config_section = malloc(sizeof(*config_section));
	config_key = malloc(sizeof(*config_key));
	config_params = malloc(sizeof(*config_params));
	config_section[config_num] = malloc(1024 * sizeof(char));
	config_key[config_num] = malloc(1024 * sizeof(char));
	config_params[config_num] = malloc(1024 * sizeof(char));
	
    while ((lread = getline(&line, &len, fp)) != -1) {
        if (!strncmp(line,"[",1)){
        	strcpy(temp_section,line);
			temp_section++;
			temp_section[strlen(temp_section)-2]=0;
        }
		if (strchr(line,'=') && strncmp(line,"#",1)){
			char *token = strtok(line, " =\n");
			temp_key = token;
			token = strtok(NULL, " =\n");
			temp_params = token;
			strcpy(config_section[config_num],temp_section);
			strcpy(config_key[config_num],temp_key);
			strcpy(config_params[config_num],temp_params);
			config_num++;
			config_section = realloc(config_section, (config_num+1)*sizeof(*config_section));
			config_key = realloc(config_key, (config_num+1)*sizeof(*config_key));
			config_params = realloc(config_params, (config_num+1)*sizeof(*config_params));
			config_section[config_num] = malloc(1024 * sizeof(char));
			config_key[config_num] = malloc(1024 * sizeof(char));
			config_params[config_num] = malloc(1024 * sizeof(char));
		}
    }

}