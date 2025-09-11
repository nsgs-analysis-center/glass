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

#include "glass_utils.h"

void ptmcmc(struct Model **model, struct Chain *chain, struct Flags *flags)
{
    int a, b;
    int olda, oldb;
    
    double heat1, heat2;
    double logL1, logL2;
    double dlogL;
    double H;
    double alpha;
    double beta;
    
    int NC = chain->NC;
    
    //b = (int)(ran2(seed)*((double)(chain->NC-1)));
    for(b=NC-1; b>0; b--)
    {
        a = b - 1;
        chain->acceptance[a]=0;
        
        olda = chain->index[a];
        oldb = chain->index[b];
        
        heat1 = chain->temperature[a];
        heat2 = chain->temperature[b];
        
        logL1 = model[olda]->logL + model[olda]->logLnorm;
        logL2 = model[oldb]->logL + model[oldb]->logLnorm;
        
        //Hot chains jump more rarely
        if(rand_r_U_0_1(&chain->r[a])<1.0)
        {
            dlogL = logL2 - logL1;
            H  = (heat2 - heat1)/(heat2*heat1);
            
            alpha = exp(dlogL*H);
            beta  = rand_r_U_0_1(&chain->r[a]);
            
            if(alpha >= beta)
            {
                chain->index[a] = oldb;
                chain->index[b] = olda;
                chain->acceptance[a]=1;
            }
        }
    }
}

void adapt_temperature_ladder(struct Chain *chain, int mcmc)
{
    int ic;
    
    int NC = chain->NC;
    
    double S[NC];
    double A[NC][2];
    
    double nu=10;
    double t0=10000.;
    
    for(ic=1; ic<NC-1; ic++)
    {
        S[ic] = log(chain->temperature[ic] - chain->temperature[ic-1]);
        A[ic][0] = chain->acceptance[ic-1];
        A[ic][1] = chain->acceptance[ic];
    }
    
    for(ic=1; ic<NC-1; ic++)
    {
        S[ic] += (A[ic][0] - A[ic][1])*(t0/((double)mcmc+t0))/nu;
        
        chain->temperature[ic] = chain->temperature[ic-1] + exp(S[ic]);
        
        if(chain->temperature[ic]/chain->temperature[ic-1] < 1.1) chain->temperature[ic] = chain->temperature[ic-1]*1.1;
    }//end loop over ic
}//end adapt function

double draw_from_uniform_prior(UNUSED struct Data *data, struct Model *model, UNUSED struct Source *source, UNUSED struct Proposal *proposal, double *params, unsigned int *seed)
{
    double logQ = 0.0;
    
    for(int n=0; n<model->NP; n++)
    {
        params[n] = model->prior[n][0] + rand_r_U_0_1(seed)*(model->prior[n][1]-model->prior[n][0]);
        logQ -= model->logPriorVolume[n];
    }
    
    return logQ;
}

double uniform_prior_density(struct Data *data, struct Model *model, UNUSED struct Source *source, UNUSED struct Proposal *proposal, double *params)
{
    double logQ = 0.0;
    
    for(int n=0; n<model->NP; n++) logQ -= model->logPriorVolume[n];
    
    return logQ;
}

double symmetric_density(UNUSED struct Data *data, UNUSED struct Model *model, UNUSED struct Source *source, UNUSED struct Proposal *proposal, UNUSED double *params)
{
    return 0.0;
}

void print_acceptance_rates(struct Proposal **proposal, int NProp, int ic, FILE *fptr)
{
    fprintf(fptr,"Acceptance rates for chain %i:\n", ic);
    fprintf(fptr," MCMC\n");
    for(int n=0; n<NProp; n++)
    {
        if(proposal[n]->weight > 0) fprintf(fptr,"   %.1e  [%s]\n", (double)proposal[n]->accept[ic]/(double)proposal[n]->trial[ic],proposal[n]->name);
    }
    fprintf(fptr," RJMCMC\n");
    for(int n=0; n<NProp; n++)
    {
        if(proposal[n]->rjweight > 0) fprintf(fptr,"   %.1e  [%s]\n", (double)proposal[n]->accept[ic]/(double)proposal[n]->trial[ic],proposal[n]->name);
    }
}
