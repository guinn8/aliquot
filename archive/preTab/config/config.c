#include "config.h"
const char configPath[] = "../data/config.cfg";
void readConfig(int count, int * ret){
    config_t cfg;
    
    

    config_init(&cfg);



    /* Read the file. If there is an error, report it and exit. */
    if(! config_read_file(&cfg, configPath)){
        fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg), config_error_line(&cfg), config_error_text(&cfg));
        config_destroy(&cfg);
      
    }
    
    //get args and catch errors
    if(!config_lookup_int(&cfg, "blocksize", &ret[0])) fprintf(stderr, "No 'blocksize' setting in configuration file.\n");
    if(!config_lookup_int(&cfg, "files_start", &ret[1])) fprintf(stderr, "No 'files_start' setting in configuration file.\n");
    if(!config_lookup_int(&cfg, "files_end", &ret[2])) fprintf(stderr, "No 'files_end' setting in configuration file.\n");
    if(!config_lookup_int(&cfg, "sigma_files", &ret[3])) fprintf(stderr, "No 'sigma_files' setting in configuration file.\n");
    if(!config_lookup_int(&cfg, "max_solutions", &ret[4])) fprintf(stderr, "No 'max_solutions' setting in configuration file.\n");
    if(!config_lookup_int(&cfg, "bound_0", &ret[5])) fprintf(stderr, "No 'bound_0' setting in configuration file.\n");
            
    config_destroy(&cfg);
   
}

int writeConfig(const long unsigned args[]){

	static const char *output_file = configPath;
  	config_t cfg;
	config_setting_t *root, *setting;
	
	
	config_init(&cfg);
	root = config_root_setting(&cfg);

	setting = config_setting_add(root, "blocksize", CONFIG_TYPE_INT);
	config_setting_set_int(setting, args[0]);

	setting = config_setting_add(root, "files_start", CONFIG_TYPE_INT);
	config_setting_set_int(setting, args[1]);

	setting = config_setting_add(root, "files_end", CONFIG_TYPE_INT);
	config_setting_set_int(setting, args[2]);

	setting = config_setting_add(root, "sigma_files", CONFIG_TYPE_INT);
	config_setting_set_int(setting, args[3]);

	setting = config_setting_add(root, "max_solutions", CONFIG_TYPE_INT);
	config_setting_set_int(setting, args[4]);

	setting = config_setting_add(root, "bound_0", CONFIG_TYPE_INT);
	config_setting_set_int(setting, args[5]);

	/* Write out the new configuration. */
	if(! config_write_file(&cfg, output_file))
	{
		fprintf(stderr, "Error while writing file.\n");
		config_destroy(&cfg);
		return(EXIT_FAILURE);
	}

	fprintf(stderr, "New configuration successfully written to: %s\n",
			output_file);

	config_destroy(&cfg);
	return(EXIT_SUCCESS);
}



