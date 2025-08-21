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
#include <glass_noise.h>
#include "glass_mbh.h"

void mbh_mcmc(struct Orbit *orbit, struct Data *data, struct Model *model, struct Model *trial, struct Chain *chain, struct Flags *flags, struct Prior *prior, struct Proposal **proposal, int ic)
{
    double logH  = 0.0; //(log) Hastings ratio
    double loga  = 1.0; //(log) transition probability
    
    double logPx  = 0.0; //(log) prior density for model x (current state)
    double logPy  = 0.0; //(log) prior density for model y (proposed state)
    double logQyx = 0.0; //(log) proposal denstiy from x->y
    double logQxy = 0.0; //(log) proposal density from y->x
    
    //shorthand pointers
    struct Model *model_x = model;
    struct Model *model_y = trial;
    
    copy_model_lite(model_x,model_y);
    
    //pick a source to update
    int n = (int)(rand_r_U_0_1(&chain->r[ic])*(double)model_x->Nlive);
        
    //more shorthand pointers
    struct Source *source_x = model_x->source[n];
    struct Source *source_y = model_y->source[n];

    //choose proposal distribution
    double draw;
    int nprop;
    do
    {
        nprop = (int)floor((MBH_PROPOSAL_NPROP)*rand_r_U_0_1(&chain->r[ic]));
        draw = rand_r_U_0_1(&chain->r[ic]);
    }while(proposal[nprop]->weight <= draw);
    
    proposal[nprop]->trial[ic]++;

    //call proposal function to update source parameters
    (*proposal[nprop]->function)(data, model_x, source_y, proposal[nprop], source_y->params, &chain->r[ic]);

    //call associated proposal density functions
    logQyx = (*proposal[nprop]->density)(data, model_x, source_y, proposal[nprop], source_y->params);
    logQxy = (*proposal[nprop]->density)(data, model_x, source_x, proposal[nprop], source_x->params);

    //get priors for x and y
    logPx = (*prior->density)(flags, data, model_x, prior, source_x->params);
    logPy = (*prior->density)(flags, data, model_y, prior, source_y->params);

    if(logPy > -INFINITY)
    {
        if(!flags->prior)
        {
            //form master template
            (*model_y->generate_signal_model)(orbit, data, model_y, 0);
            
            //get likelihood for y using delta log likelihood method
            model_y->logL  = gaussian_log_likelhood_wavelet(data, model_y);

            /*compare fast and slow likelihood
            double dlogLx = delta_log_likelihood(data, model_x);
            double dlogLy = delta_log_likelihood(data, model_y);
            printf("logLx=%.12g, logLy=%.12g, dlogL=%.12g, dlogLx=%.12g, dlogLy=%.12g, dlogLtemp=%.12g\n",
                   model_x->logL,model_y->logL,
                   model_y->logL-model_x->logL,
                   dlogLx,dlogLy,
                   dlogLy-dlogLx);*/
            
            /*
             H = [p(d|y)/p(d|x)]/T x p(y)/p(x) x q(x|y)/q(y|x)
             */
            logH += (model_y->logL - model_x->logL)/chain->temperature[ic]; //delta logL
        }
        logH += logPy  - logPx;  //priors
        logH += logQxy - logQyx; //proposals
        
        loga = log(rand_r_U_0_1(&chain->r[ic]));
        
        if(isfinite(logH) && logH > loga)
        {

            proposal[nprop]->accept[ic]++;
            copy_model_lite(model_y,model_x);

        }
    }
}

void initialize_mbh_state(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, struct Proposal **proposal, struct Model **model, struct Model **trial, struct Source **inj_vec)
{
    fprintf(stdout,"\n======== Initialize MBH Sampler State =======\n");

    int NC = chain->NC;
    int DMAX = flags->DMAX;
    
    //keep track of number of injections
    struct Source *inj;
    int Ninj=0;
    
    for(int ic=0; ic<NC; ic++)
    {
        trial[ic] = malloc(sizeof(struct Model));
        model[ic] = malloc(sizeof(struct Model));
        
        alloc_model(data,trial[ic],MBH_MODEL_NP,DMAX);
        alloc_model(data,model[ic],MBH_MODEL_NP,DMAX);
        
        model[ic]->log_likelihood = &delta_log_likelihood;
        trial[ic]->log_likelihood = &delta_log_likelihood;
        
        model[ic]->generate_signal_model = &generate_mbh_signal_model;
        trial[ic]->generate_signal_model = &generate_mbh_signal_model;
        
        if(ic==0)set_uniform_mbh_prior(flags, model[ic], data, 1);
        else     set_uniform_mbh_prior(flags, model[ic], data, 0);
                
        //override noise model w/ stationary version
        if(flags->stationary) GetStationaryNoiseModel(data, orbit, flags, data->noise);

        //set noise model
        copy_noise(data->noise, model[ic]->noise);
        
        //draw signal model
        for(int n=0; n<DMAX; n++)
        {
            if(flags->cheat && inj_vec[n]->tref>data->t0)
            {
                
                if(ic==0)Ninj++;
                inj=inj_vec[n];
                
                //map parameters to vector
                model[ic]->source[n]->Mchirp   = inj->Mchirp;
                model[ic]->source[n]->Mtotal   = inj->Mtotal;
                model[ic]->source[n]->chi1     = inj->chi1;
                model[ic]->source[n]->chi2     = inj->chi2;
                model[ic]->source[n]->phiref   = inj->phiref;
                model[ic]->source[n]->tref     = inj->tref;
                model[ic]->source[n]->dL       = inj->dL;
                model[ic]->source[n]->costheta = inj->costheta;
                model[ic]->source[n]->phi      = inj->phi;
                model[ic]->source[n]->psi      = inj->psi;
                model[ic]->source[n]->cosi     = inj->cosi;
                map_mbh_params_to_array(model[ic]->source[n], model[ic]->source[n]->params);
                if(ic==0 && n==0 && flags->stationary)
                    fprintf(stdout,"   ...stationary SNR for source %i=%g\n",n,snr_wavelet(inj,model[ic]->noise));
            }
            else
            {
                do {
                    draw_from_mbh_prior(data, model[ic], model[ic]->source[n], proposal[0], model[ic]->source[n]->params, &chain->r[ic]);
                } while (check_mbh_prior(model[ic]->source[n]->params,model[ic]->prior));
                
                double m1 = 9.840620473718e+05;
                double m2 = 8.212835097523e+05;
                
                model[ic]->source[n]->Mchirp   = chirpmass(m1,m2);
                model[ic]->source[n]->Mtotal   = m1+m2;
                model[ic]->source[n]->chi1     = 4.525488035553191e-01;
                model[ic]->source[n]->chi2     = 5.581335050668826e-01;
                model[ic]->source[n]->phiref   = 2.982619790826792e+00;
                model[ic]->source[n]->tref     = 4.800035465863327e+06;
                model[ic]->source[n]->dL       = 1.763565935270e+01;
                model[ic]->source[n]->costheta = -4.741945333726355e-01;
                model[ic]->source[n]->phi      = 5.823731551490122e-01;
                model[ic]->source[n]->psi      = 1.840857645136265e-01;
                model[ic]->source[n]->cosi     = -2.544566370040902e-02;
                map_mbh_params_to_array(model[ic]->source[n], model[ic]->source[n]->params);
                
            }
            map_array_to_mbh_params(model[ic]->source[n], model[ic]->source[n]->params);

            mbh_fisher(orbit, data, model[ic]->source[n], data->noise);

            model[ic]->source[n]->fisher_update_flag=0;
        }
        
        //initialize sampler to proper size of model
        if(flags->cheat)
        {
            model[ic]->Neff = Ninj+1;
            model[ic]->Nlive= Ninj;
        }
        
        /*
         Form master model & compute likelihood of starting position
         */
        
        //signal
        (*model[ic]->generate_signal_model)(orbit,data,model[ic],-1);

        //noise
        for(int i=0; i<data->Nchannel; i++)
        {
            for(int j=i; j<data->Nchannel; j++)
            {
                memcpy(model[ic]->noise->C[i][j], data->noise->C[i][j], data->N*sizeof(double));
                memcpy(model[ic]->noise->C[j][i], data->noise->C[i][j], data->N*sizeof(double));
            }
        }

        invert_noise_covariance_matrix(model[ic]->noise);
        
        
        /*calibration error
        if(flags->calibration)
        {
            draw_calibration_parameters(data, model[ic], &chain->r[ic]);
            generate_calibration_model(data, model[ic]);
            apply_calibration_model(data, model[ic]);
        }*/
        
        //likelihood
        if(!flags->prior)
        {
            model[ic]->logL = gaussian_log_likelihood(data, model[ic]);
            model[ic]->logLnorm = gaussian_log_likelihood_model_norm(data,model[ic]);
        }
        else model[ic]->logL = model[ic]->logLnorm = 0.0;
        
        printf("  log likelihood of chain %i: %g\n",ic,model[ic]->logL);
        
        if(ic==0) chain->logLmax += model[ic]->logL + model[ic]->logLnorm;
        
    }//end loop over chains
    
    /* initialize differential evolution buffer */
    //first remember which proposal is for differential evoltuion
    int nprop = -1;
    do {
        nprop++;
    } while (strcmp(proposal[nprop]->name,"de jump"));
    
    //fill buffer with prior draws
    for(int n=0; n<proposal[nprop]->size; n++)
    {
        double *params = proposal[nprop]->matrix[n];
        do {
            draw_from_mbh_prior(data,model[0],model[0]->source[0],proposal[nprop],params,&chain->r[0]);
        } while (check_mbh_prior(params,model[0]->prior));
        
    }

    fprintf(stdout,"=============================================\n");
}
