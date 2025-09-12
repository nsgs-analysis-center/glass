/*
 * Copyright 2025 Tyson B. Littenberg & Neil Cornish
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

double draw_from_mbh_prior(struct Data *data, struct Model *model, struct Source * source, struct Proposal *proposal, double *params, unsigned int *seed)
{
    return draw_from_uniform_prior(data,model,source,proposal,params,seed);
}

void setup_differential_evolution_proposal(struct Proposal *proposal)
{
    //differential evolution buffer
    proposal->size = MBH_DIFFERENTIAL_EVOLUTION_BUFFER;
    proposal->matrix = double_matrix(proposal->size, MBH_MODEL_NP);
}

double differential_evolution_jump(struct Data *data, struct Model *model, struct Source * source, struct Proposal *proposal, double *params, unsigned int *seed)
{
    // aliases to contents of the proposal structure
    double **history = proposal->matrix;
    double N_history = proposal->size;

    // pick two points from the history
    int i = (int)(N_history * rand_r_U_0_1(seed));
    int j = (int)(N_history * rand_r_U_0_1(seed));
    while(i==j) j = (int)(N_history * rand_r_U_0_1(seed));

    // scale of jump
    double alpha = 1.0;
    if(rand_r_U_0_1(seed)<0.9) alpha = rand_r_N_0_1(seed)*0.5;

    // jump in all parameters
    if(rand_r_U_0_1(seed)<0.5)
    {
        for(int n=0; n<MBH_MODEL_NP; n++) 
            params[n] += alpha*(history[i][n] - history[j][n]);
    }
    // just intrinsic parameters
    else
    {
        params[0] = params[0]+alpha*(history[i][0]-history[j][0]);
        params[1] = params[1]+alpha*(history[i][1]-history[j][1]);
        params[2] = params[2]+alpha*(history[i][2]-history[j][2]);
        params[3] = params[3]+alpha*(history[i][3]-history[j][3]);
    }

    //differential evolution jump is symmetric
    return 0.0;
}

void update_differential_evolution_buffer(struct Proposal **proposal, struct Model *model, unsigned int *seed)
{
        
    //first remember which proposal is for differential evoltuion
    int nprop = -1;
    do {
        nprop++;
    } while (strcmp(proposal[nprop]->name,"de jump"));
    
    double **buffer = proposal[nprop]->matrix;
    int Nbuffer = proposal[nprop]->size;
    
    //source parameters from current state of model
    double *params = model->source[0]->params;
    
    //pick a random sample from the buffer to replace
    int n = (int)(Nbuffer * rand_r_U_0_1(seed));
    
    memcpy(buffer[n],params,MBH_MODEL_NP*sizeof(double));
}

double draw_from_mbh_fisher(UNUSED struct Data *data, struct Model *model, struct Source *source, UNUSED struct Proposal *proposal, double *params, unsigned int *seed)
{
    int i,j;
    //double sqNP = sqrt((double)source->NP);
    double Amps[MBH_MODEL_NP];
    double jump[MBH_MODEL_NP];
    
    //draw the eigen-jump amplitudes from N[0,1] scaled by evalue & dimension
    for(i=0; i<MBH_MODEL_NP; i++)
    {
        Amps[i] = rand_r_N_0_1(seed)/sqrt(source->fisher_evalue[i]);
        jump[i] = 0.0;
    }
    
    //choose one eigenvector to jump along
    i = (int)(rand_r_U_0_1(seed)*(double)MBH_MODEL_NP);
    for (j=0; j<MBH_MODEL_NP; j++)
    {
        jump[j] += Amps[i]*source->fisher_evectr[j][i];
    }

    //check jump value, set to small value if singular
    for(i=0; i<MBH_MODEL_NP; i++) if(!isfinite(jump[i])) jump[i] = 0.01*source->params[i];
    
    //jump from current position
    for(i=0; i<MBH_MODEL_NP; i++) params[i] = source->params[i] + jump[i];
    
    //safety check for cos(latitude) parameters
    //cosine co-latitude
    if(params[7] >=  1.) params[7] = source->params[7] - jump[7];
    if(params[7] <= -1.) params[7] = source->params[7] - jump[7];
    //cosine inclination
    if(params[10] >=  1.) params[10] = source->params[10] - jump[10];
    if(params[10] <= -1.) params[10] = source->params[10] - jump[10];

    for(int j=0; j<MBH_MODEL_NP; j++)
    {
        if(params[j]!=params[j])
        {
            fprintf(stderr,"draw_from_fisher: params[%i]=%g, N[%g,%g]\n",j,params[j],source->params[j],jump[j]);
            return -INFINITY;
        }
    }
        
    //not updating Fisher between moves, proposal is symmetric
    return 0.0;
}

void initialize_mbh_proposal(struct Orbit *orbit, struct Data *data, struct Prior *prior, struct Chain *chain, struct Flags *flags, struct Catalog *catalog, struct Proposal **proposal, int NMAX)
{
    fprintf(stdout,"\n===== Initialize MBH Proposal Structure =====\n");

    int NC = chain->NC;
    double check  =0.0;
    double rjcheck=0.0;
    
    for(int i=0; i<MBH_PROPOSAL_NPROP; i++)
    {
        proposal[i] = malloc(sizeof(struct Proposal));

        proposal[i]->trial  = malloc(NC*sizeof(int));
        proposal[i]->accept = malloc(NC*sizeof(int));
        
        for(int ic=0; ic<NC; ic++)
        {
            proposal[i]->trial[ic]  = 1;
            proposal[i]->accept[ic] = 0;
        }
        
        switch(i)
        {
            case 0:
                sprintf(proposal[i]->name,"prior");
                proposal[i]->function = &draw_from_mbh_prior;
                proposal[i]->density  = &uniform_prior_density;
                proposal[i]->weight   = 0.2;
                proposal[i]->rjweight = 1.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 1:
                sprintf(proposal[i]->name,"de jump");
                setup_differential_evolution_proposal(proposal[i]);
                proposal[i]->function = &differential_evolution_jump;
                proposal[i]->density  = &symmetric_density;
                proposal[i]->weight   = 0.4;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
                break;
            case 2:
                sprintf(proposal[i]->name, "fisher");
                proposal[i]->function = &draw_from_mbh_fisher;
                proposal[i]->density  = &symmetric_density;
                proposal[i]->weight   = 0.4;
                proposal[i]->rjweight = 0.0;
                check   += proposal[i]->weight;
                rjcheck += proposal[i]->rjweight;
            default:
                break;
        }
    }
    
    if(check!=1.0 || rjcheck!=1.0)
    {
        printf("Error: Set of proposal distributions is not normalized\n");
        printf("       check   = %g\n",check);
        printf("       rjcheck = %g\n",check);
        exit(1);
    }
    fprintf(stdout,"=============================================\n");
}
