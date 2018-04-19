#ifndef __CONFIGPARSER__
#define __CONFIGPARSER__

int config_num;
char ** config_section;
char ** config_key;
char ** config_params;

#ifdef __cplusplus
extern "C" void read_config(char *filename);
extern "C" int config_get_int(char *search_section, char *search_key);
extern "C" float config_get_float(char *search_section, char *search_key);
extern "C" char *config_get_string(char *search_section, char *search_key);
extern "C" int *config_get_intarr(char *search_section, char *search_key, int *size);
#endif

#endif
