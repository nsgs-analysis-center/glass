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

}
void parse_mbh_args(int argc, char **argv, struct Flags *flags)
{
    //copy argv since getopt permutes order
    char **argv_copy=malloc((argc+1) * sizeof *argv_copy);
    copy_argv(argc,argv,argv_copy);
    opterr=0; //suppress warnings about unknown arguments

    flags->DMAX = 1;

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
    print_mbh_chain_state(data, chain, model[n], flags, chain->chainFile[0], step);
    
    
    //Print sampling parameters
    int D = model[n]->Nlive;
    for(i=0; i<D; i++)
    {
        print_mbh_source_params(data,model[n]->source[i],chain->parameterFile[0]);
        fprintf(chain->parameterFile[0],"\n");
        print_mbh_source_params(data,model[n]->source[i],stdout);
        fprintf(stdout,"\n");
        fflush(chain->parameterFile[0]);
        fflush(stdout);
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

void print_mbh_chain_state(struct Data *data, struct Chain *chain, struct Model *model, struct Flags *flags, FILE *fptr, int step)
{
    fprintf(fptr, "%i ",step);
    fprintf(fptr, "%i ",model->Nlive);
    fprintf(fptr, "%lg ",model->logL);
    fprintf(fptr, "\n");
}
void scan_mbh_chain_state(struct Data *data, struct Chain *chain, struct Model *model, struct Flags *flags, FILE *fptr, int *step)
{
    int check = 0;
    check += fscanf(fptr, "%i",step);
    check += fscanf(fptr, "%i",&model->Nlive);
    check += fscanf(fptr, "%lg",&model->logL);

    if(!check)
    {
        fprintf(stderr,"Error reading checkpoint files\n");
        exit(1);
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
