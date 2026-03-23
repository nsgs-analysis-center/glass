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

/**
 \brief Converts array expected by glass_mbh_waveform.c to
 physical MBH parameters
 */
void map_array_to_mbh_params(struct Source *source, double *params)
{
    source->Mchirp   = exp(params[0]);
    source->Mtotal   = exp(params[1]);
    source->chi1     = params[2];
    source->chi2     = params[3];
    source->phiref   = params[4];
    source->tref     = params[5];
    source->dL       = exp(params[6]);
    source->costheta = params[7];
    source->phi      = params[8];
    source->psi      = params[9];
    source->cosi     = params[10];
}

/**
 \brief Converts physical MBH parameters to array expected by glass_mbh_waveform.c
 */
void map_mbh_params_to_array(struct Source *source, double *params)
{
    params[0]  = log(source->Mchirp);
    params[1]  = log(source->Mtotal);
    params[2]  = source->chi1;
    params[3]  = source->chi2;
    params[4]  = source->phiref;
    params[5]  = source->tref;
    params[6]  = log(source->dL);
    params[7]  = source->costheta;
    params[8]  = source->phi;
    params[9]  = source->psi;
    params[10] = source->cosi;
}

void generate_mbh_signal_model(struct Orbit *orbit, struct Data *data, struct Model *model, int source_id)
{
    int i,n;

    /*
     Generate waveform for Source
     */
    
    struct Source *source = model->source[source_id];
    
    for(i=0; i<data->N; i++)
    {
        source->tdi->X[i]=0.0;
        source->tdi->Y[i]=0.0;
        source->tdi->Z[i]=0.0;
        source->list[i]=0;
    }
    
    map_array_to_mbh_params(source, source->params);
    
    //Simulate gravitational wave signal
    mbh_fd_waveform(orbit, data->wdm, data->T, model->t0, source->params, source->list, &source->Nlist, source->tdi->X, source->tdi->Y, source->tdi->Z);
    
    
    /*
     Add Source to Model
     */
    for(n=0; n<data->N; n++)
    {
        model->tdi->X[n]=0.0;
        model->tdi->Y[n]=0.0;
        model->tdi->Z[n]=0.0;
        model->list[n]=0;
    }
    model->Nlist = source->Nlist;

    for(int n=0; n<source->Nlist; n++)
    {
        int k = source->list[n];
        if(k>=0 && k<data->N)
        {
            model->tdi->X[k] += source->tdi->X[k];
            model->tdi->Y[k] += source->tdi->Y[k];
            model->tdi->Z[k] += source->tdi->Z[k];
            model->list[n] = k;
        }
    }

    //get union of list
    //this is superfluous unless we end up with multiple live MBH's per model structure
    //list_union(model->list, source->list, model->Nlist, source->Nlist, model->list, &model->Nlist);
}
