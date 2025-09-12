/*
 * Copyright 2025 Tyson B. Littenberg & Neil J. Cornish
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

int mbh_mass_ratio_check(double m1, double m2)
{
    double q = m2/m1;
    if(q>MBH_MAX_MASS_RATIO || q<MBH_MIN_MASS_RATIO) return 1;
    else return 0;
}

void set_uniform_mbh_prior(struct Flags *flags, struct Model *model, struct Data *data, int verbose)
{
    
    model->t0 = data->t0; //TODO move this to model setup
    
    //TODO: assign priors by parameter name, use mapper to get into vector (more robust to changes)
    
    //log chirp mass
    model->prior[0][0] = log(1.0e2);
    model->prior[0][1] = log(0.44*5.0e8);
    
    //log total mass
    model->prior[1][0] = log(1.0e3);
    model->prior[1][1] = log(5.0e8);
    
    //chi1
    model->prior[2][0] = -0.999;
    model->prior[2][1] = +0.999;
    
    //chi2
    model->prior[3][0] = -0.999;
    model->prior[3][1] = +0.999;
    
    //merger phase
    model->prior[4][0] = 0.0;
    model->prior[4][1] = PI2;

    //merger time
    if(flags->triggerTime < 0.0)
    {
        model->prior[5][0] = 1.01*data->t0; //TODO: 1.01 seems pretty arbitrary
        model->prior[5][1] = 2*data->T;     //TODO: max tc is too restrictive
        model->prior[5][1] = 0.99*data->T;     //TODO: max tc is too restrictive
    }
    else
    {
        model->prior[5][0] = flags->triggerTime - DAY;
        model->prior[5][1] = flags->triggerTime + DAY;
    }
    //log distance (Gpc)
    model->prior[6][0] = log(0.1);
    model->prior[6][1] = log(1.0e3);

    //cos EclipticCoLatitude
    model->prior[7][0] = -1.0;
    model->prior[7][1] =  1.0;

    //EclipticLongitude
    model->prior[8][0] = 0.0;
    model->prior[8][1] = PI2;

    //polarization angle
    model->prior[9][0] = 0.0;
    model->prior[9][1] = M_PI;

    //cos inclination
    model->prior[10][0] = -1.0;
    model->prior[10][1] =  1.0;

        
    if(verbose && !flags->quiet)
    {
        fprintf(stdout,"  MBH priors:\n");
    }
    
    
    //set prior volume
    for(int n=0; n<MBH_MODEL_NP; n++) model->logPriorVolume[n] = log(model->prior[n][1]-model->prior[n][0]);
    
}

int check_mbh_prior(double *params, double **uniform_prior)
{
    //nan check
    for(int n=0; n<MBH_MODEL_NP; n++) if(params[n]!=params[n]) return 1;
    
    //Chirp mass
    if(params[0]<uniform_prior[0][0] || params[0]>uniform_prior[0][1]) return 1;
    
    //Total mass
    if(params[1]<uniform_prior[1][0] || params[1]>uniform_prior[1][1]) return 1;

    //Spin 1
    if(params[2]<uniform_prior[2][0] || params[2]>uniform_prior[2][1]) return 1;
    
    //Spin 2
    if(params[3]<uniform_prior[3][0] || params[3]>uniform_prior[3][1]) return 1;

    //Merger phase
    while(params[4]<uniform_prior[4][0]) params[4] += PI2;
    while(params[4]>uniform_prior[4][1]) params[4] -= PI2;

    //Merger time
    if(params[5]<uniform_prior[5][0] || params[5]>uniform_prior[5][1]) return 1;

    //Distance
    if(params[6]<uniform_prior[6][0] || params[6]>uniform_prior[6][1]) return 1;

    //Cosine co-latitude
    if(params[7]<uniform_prior[7][0] || params[7]>uniform_prior[7][1]) return 1;

    //longitude
    while(params[8]<uniform_prior[8][0]) params[8] += PI2;
    while(params[8]>uniform_prior[8][1]) params[8] -= PI2;

    //polarization
    while(params[9]<uniform_prior[9][0]) params[9] += M_PI;
    while(params[9]>uniform_prior[9][1]) params[9] -= M_PI;

    //cosine inclination
    if(params[10]<uniform_prior[10][0] || params[10]>uniform_prior[10][1]) return 1;
    
    //extrajudicial check on component masses
    double Mchirp,Mtotal,m1,m2;
    Mchirp = exp(params[0]);
    Mtotal = exp(params[1]);
    component_masses(Mchirp,Mtotal,&m1,&m2);
    if(mbh_mass_ratio_check(m1,m2)) return 1;

    return 0;
}

double evaluate_mbh_prior(struct Flags *flags, struct Data *data, struct Model *model, struct Prior *prior, double *params)
{
    double logP=0.0;
    
    // check that everything is in range
    if(check_mbh_prior(params, model->prior)) return -INFINITY;
    
    // if so get prior volume
    for(int n=0; n<MBH_MODEL_NP; n++) logP -= model->logPriorVolume[n];
    
    return logP;
}
