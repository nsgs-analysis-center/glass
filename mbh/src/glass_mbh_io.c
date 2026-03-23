/*
 * Copyright 2025 Tyson B. Littenberg
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <glass_utils.h>
#include "glass_mbh.h"

void print_mbh_usage()
{
    fprintf(stdout,"\n");
    fprintf(stdout,"================ MBH Usage: =============== \n");
    fprintf(stdout,"REQUIRED:\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"OPTIONAL:\n");

    //Priors & Proposals
    fprintf(stdout,"       ==== Priors & Proposals ==== \n");
    fprintf(stdout,"       --trigger     : event trigger time for analysis     \n");

    //Injections
    fprintf(stdout,"       ======== Injections ======== \n");
    fprintf(stdout,"       --inj         : inject signal                       \n");
    fprintf(stdout,"       --injseed     : seed for injection parameters       \n");
    fprintf(stdout,"\n");

}
void parse_mbh_args(int argc, char **argv, struct Flags *flags)
{
    //copy argv since getopt permutes order
    char **argv_copy=malloc((argc+1) * sizeof *argv_copy);
    copy_argv(argc,argv,argv_copy);
    opterr=0; //suppress warnings about unknown arguments

    //Set defaults
    flags->NINJ  = 0;
    flags->DMAX  = 1;
    flags->cheat = 0;
    flags->triggerTime = -1.0;

    flags->injFile = malloc(flags->DMAX *sizeof(char *));
    for(int n=0; n<flags->DMAX ; n++) flags->injFile[n] = malloc(1024*sizeof(char));
    
    //Specifying the expected options
    static struct option long_options[] =
    {
        /* These options set a flag. */
        {"injseed",    required_argument, 0, 0},
        {"inj",        required_argument, 0, 0},
        {"trigger",    required_argument, 0, 0},
        
        /* These options don’t set a flag.
         We distinguish them by their indices. */
        {"cheat",       no_argument, 0, 0 },
        {0, 0, 0, 0}
    };
    
    
    int opt=0;
    int long_index=0;
    
    //Loop through argv string and pluck out arguments
    while ((opt = getopt_long_only(argc, argv_copy,"apl:b:", long_options, &long_index )) != -1)
    {
        switch (opt)
        {
                
            case 0:
                if(strcmp("cheat", long_options[long_index].name) == 0) flags->cheat = 1;
                if(strcmp("trigger",long_options[long_index].name) == 0) flags->triggerTime = atof(optarg);
                if(strcmp("inj", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    sprintf(flags->injFile[flags->NINJ],"%s",optarg);
                    flags->NINJ++;
                    if(flags->NINJ>flags->DMAX )
                    {
                        fprintf(stderr,"WARNING: Requested number of injections is too large (%i/%i)\n",flags->NINJ,flags->DMAX );
                        fprintf(stderr,"Should you remove at least %i --inj arguments?\n",flags->NINJ-flags->DMAX );
                        fprintf(stderr,"Now exiting to system\n");
                        exit(1);
                    }
                }
                break;
            case 'h' :
                print_mbh_usage();
                break;
            default:
                break;
        }
    }

    
    
    //reset opt counter
    optind = 0;

    //free placeholder for argvs
    for(int i=0; i<=argc; i++)free(argv_copy[i]);
    free(argv_copy);
}

void print_mbh_chain_files(struct Data *data, struct Model **model, struct Chain *chain, struct Flags *flags, int step)
{
    int i,n,ic;
    
    //Print logL & temperature chains
    if(!flags->quiet)
    {
        fprintf(chain->likelihoodFile,  "%i ",step);
        fprintf(chain->temperatureFile, "%i ",step);
        double logL;
        for(ic=0; ic<chain->NC; ic++)
        {
            n = chain->index[ic];
            logL=0.0;
            logL += model[n]->logL+model[n]->logLnorm;
            fprintf(chain->likelihoodFile,  "%lg ",logL);
            fprintf(chain->temperatureFile, "%lg ",1./chain->temperature[ic]);
        }
        fprintf(chain->likelihoodFile, "\n");
        fprintf(chain->temperatureFile,"\n");
    }
    
    //Print cold chains
    n = chain->index[0];
    print_chain_state(data, chain, model[n], flags, chain->chainFile[0], step);
    
    
    //Print sampling parameters
    int D = model[n]->Nlive;
    for(i=0; i<D; i++)
    {
        print_mbh_source_params(data,model[n]->source[i],chain->parameterFile[0]);
        fprintf(chain->parameterFile[0],"\n");
        
//        fprintf(stdout,"%.12g %.12g ",model[n]->logL,sqrt(2.*(model[n]->logL-data->SNR2)));
//        print_mbh_source_params(data,model[n]->source[i],stdout);
//        fprintf(stdout,"\n");

        if(step>0)
        {
            if(chain->dimensionFile[D]==NULL)
            {
                char filename[MAXSTRINGSIZE];
                sprintf(filename,"%s/dimension_chain.dat.%i",chain->chainDir,D);
                if(flags->resume)chain->dimensionFile[D] = fopen(filename,"a");
                else             chain->dimensionFile[D] = fopen(filename,"w");
            }
            print_mbh_source_params(data,model[n]->source[i],chain->dimensionFile[D]);
            fprintf(chain->dimensionFile[D],"\n");
        }
    }
}

void print_mbh_source_params(struct Data *data, struct Source *source, FILE *fptr)
{
    //map to parameter names (just to make code readable)
    map_array_to_mbh_params(source, source->params);
    
    fprintf(fptr,"%.16g ",source->Mchirp);
    fprintf(fptr,"%.12g ",source->Mtotal);
    fprintf(fptr,"%.12g ",source->chi1);
    fprintf(fptr,"%.12g ",source->chi2);
    fprintf(fptr,"%.12g ",source->tref);
    fprintf(fptr,"%.12g ",source->phi);
    fprintf(fptr,"%.12g ",source->costheta);
    fprintf(fptr,"%.12g ",source->dL);
    fprintf(fptr,"%.12g ",source->cosi);
    fprintf(fptr,"%.12g ",source->psi);
    fprintf(fptr,"%.12g ",source->phiref);
}

void scan_mbh_source_params(struct Data *data, struct Source *source, FILE *fptr)
{
    int check = 0;

    check+=fscanf(fptr,"%lg",&source->Mchirp);
    check+=fscanf(fptr,"%lg",&source->Mtotal);
    check+=fscanf(fptr,"%lg",&source->chi1);
    check+=fscanf(fptr,"%lg",&source->chi2);
    check+=fscanf(fptr,"%lg",&source->tref);
    check+=fscanf(fptr,"%lg",&source->phi);
    check+=fscanf(fptr,"%lg",&source->costheta);
    check+=fscanf(fptr,"%lg",&source->dL);
    check+=fscanf(fptr,"%lg",&source->cosi);
    check+=fscanf(fptr,"%lg",&source->psi);
    check+=fscanf(fptr,"%lg",&source->phiref);

    
    if(!check)
    {
        fprintf(stdout,"Error reading source file\n");
        exit(1);
    }
    
    map_mbh_params_to_array(source, source->params);
}
