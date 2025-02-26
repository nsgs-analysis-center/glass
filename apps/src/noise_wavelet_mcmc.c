/*
 *  Copyright (C) 2024 Robert Rosati MSFC/UAH
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/**
 @file noise_wavelet_mcmc.c
 \brief Main function for wavelet-domain noise and SGWB sampler app `noise_wavelet_mcmc` 
 */

/*  REQUIRED LIBRARIES  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <sys/stat.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <omp.h>

#include <glass_utils.h>
#include <glass_noise.h>


static void print_usage()
{
    print_glass_usage();
    fprintf(stdout,"EXAMPLE:\n");
    fprintf(stdout,"noise_wavelet_mcmc --sim-noise --conf-noise --duration 7864320 --fmin 1e-4 --fmax 8e-3\n");
    fprintf(stdout,"\n");
    exit(0);
}

int main(int argc, char *argv[])
{
    fprintf(stdout, "\n============== NOISE WAVELET MCMC ===========\n");

    time_t start, stop;
    start = time(NULL);
    char filename[MAXSTRINGSIZE];
    
    /* check arguments */
    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    if(argc==1) print_usage();

    bool stationary_test = false;

    /* Allocate data structures */
    struct Data *data   = malloc(sizeof(struct Data));
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Chain *chain = malloc(sizeof(struct Chain));
    
    parse_data_args(argc,argv,data,orbit,flags,chain,"wavelet");
    if(flags->help) print_usage();

    /* Setup output directories for data and chain files */
    sprintf(data->dataDir,"%s/data",flags->runDir);
    sprintf(chain->chainDir,"%s/chains",flags->runDir);
    sprintf(chain->chkptDir,"%s/checkpoint",flags->runDir);
    
    mkdir(flags->runDir,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data->dataDir,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


    /* Initialize data structures */
    alloc_data(data, flags);

    /*
    Define and set up Orbit structure which contains spacecraft ephemerides
    */
    initialize_orbit(data, orbit, flags);
    initialize_interpolated_analytic_orbits(orbit, data->T, data->t0);

    /* Initialize chain structure and files */
    initialize_chain(chain, flags, &data->cseed, "a");


    // okay, for now we're in the very weird situation of not trying to fit the foreground or instrument params, **only** the SGWB params
    //
    /* Initialize Instrument Noise Model */
    printf("   ...initialize instrument noise model\n");
    struct Noise **psd = malloc(chain->NC*sizeof(struct Noise *)); // will save total scalograms here
    struct InstrumentModel **inst_model = malloc(chain->NC*sizeof(struct InstrumentModel *));
    struct InstrumentModel **inst_trial = malloc(chain->NC*sizeof(struct InstrumentModel *));
    for(int ic=0; ic<chain->NC; ic++)
    {
        psd[ic] = malloc(sizeof(struct Noise));
        alloc_noise(psd[ic], data->N, data->Nchannel); // note is data->N in wavelet not data->NFFT

        inst_model[ic] = malloc(sizeof(struct InstrumentModel));
        inst_trial[ic] = malloc(sizeof(struct InstrumentModel));
        initialize_instrument_model_wavelet(orbit, data, inst_model[ic]);
        initialize_instrument_model_wavelet(orbit, data, inst_trial[ic]);
    }
    // TODO for now this works because the noise model is stationary
    sprintf(filename,"%s/instrument_noise_model.dat",data->dataDir);
    print_noise_model(inst_model[0]->psd, filename);

    /* Initialize Galactic Foreground Model */
    // note: UCB wavelet app hides this inside initialize_ucb_state if stationary noise, and in the waveform calls otherwise.
    // See GetDynamicNoiseModel and GetStationaryNoiseModel
    // remember, plan here for now is to not sample over foreground model in wavelet basis, but still sample over noise/sgwb
    // doesn't look too hard tbh!
    if(flags->confNoise)printf("   ...initialize foreground noise model\n");
    struct ForegroundModel **conf_model = malloc(chain->NC*sizeof(struct ForegroundModel *));
    for(int ic=0; ic<chain->NC; ic++)
    {
        conf_model[ic] = malloc(sizeof(struct ForegroundModel));
        if(flags->confNoise) 
        {
           initialize_foreground_model_wavelet(orbit, data, conf_model[ic]);
        }
    }
    // TODO check if this even works
    if(flags->confNoise)
    {
        sprintf(filename,"%s/foreground_noise_model.dat",data->dataDir);
        print_noise_model(conf_model[0]->psd, filename);
    }


    /* Initialize SGWB Model */
    if(flags->sgwbTemplate>=0) printf("   ...initialize SGWB model\n");
    struct SGWBModel **sgwb_model = malloc(chain->NC*sizeof(struct SGWBModel *));
    struct SGWBModel **sgwb_trial = malloc(chain->NC*sizeof(struct SGWBModel *));
    for (int ic=0; ic<chain->NC; ic++)
    {
        sgwb_model[ic] = malloc(sizeof(struct SGWBModel));
        sgwb_trial[ic] = malloc(sizeof(struct SGWBModel));
        if(flags->sgwbTemplate>=0) 
        {
           initialize_sgwb_model(orbit, data, sgwb_model[ic], flags->sgwbTemplate);
           initialize_sgwb_model(orbit, data, sgwb_trial[ic], flags->sgwbTemplate);
        }
    }
    if(flags->sgwbTemplate>=0)
    {
        sprintf(filename,"%s/sgwb_noise_model.dat",data->dataDir);
        print_noise_model(sgwb_model[0]->psd, filename);
        // works because stationary
    }

    /* TODO: SGWB injections??? */

    /* Combine noise components to form covariance matrix */
    // TODO need stationary flag
    for(int ic=0; ic<chain->NC; ic++)
    {
        if(!flags->confNoise || flags->sgwbTemplate==0){
            printf("error: only support conf and sgwb on!");
            return -1;
        }
        generate_full_dynamic_covariance_matrix(data->wdm,inst_model[ic],conf_noise[ic],sgwb_model[ic],psd[ic]);
    }

    /* get initial likelihood */
    for(int ic=0; ic<chain->NC; ic++)
    {
        invert_noise_covariance_matrix(psds[ic]);
        psds[ic]->logL = noise_log_likelihood(data, inst_model[ic]->psd);
    }

    sprintf(filename,"%s/full_noise_model.dat",data->dataDir);
    print_noise_model(inst_model[0]->psd, filename);

    //MCMC
    printf("\n==== Noise Wavelet MCMC Sampler ====\n");
    // TODO not even slightly wavelet yet
    
    sprintf(filename,"%s/noise_chain.dat",chain->chainDir);
    FILE *noiseChainFile = fopen(filename,"w");

    FILE *foregroundChainFile = NULL;
    if(flags->confNoise)
    {
        sprintf(filename,"%s/foreground_chain.dat",chain->chainDir);
        foregroundChainFile = fopen(filename,"w");
    }
    FILE *sgwbChainFile = NULL;
    if(flags->sgwbTemplate>=0)
    {
        sprintf(filename,"%s/sgwb_chain.dat",chain->chainDir);
        sgwbChainFile = fopen(filename,"w");
    }
    
    int numThreads;
    int step = 0;
    int NC = chain->NC;
    
    #pragma omp parallel num_threads(flags->threads)
    {
        int threadID;
        
        //Save individual thread number
        threadID = omp_get_thread_num();
        
        //Only one thread runs this section
        if(threadID==0)  numThreads = omp_get_num_threads();
        
        #pragma omp barrier
        
        /* The MCMC loop */
        for(; step<flags->NMCMC;)
        {
            #pragma omp barrier
            
            // (parallel) loop over chains
            for(int ic=threadID; ic<NC; ic+=numThreads)
            {
                struct Noise *psd_ptr = psd[chain->index[ic]];
                struct InstrumentModel *inst_model_ptr = inst_model[chain->index[ic]];
                struct InstrumentModel *inst_trial_ptr = inst_trial[chain->index[ic]];
                struct ForegroundModel *conf_model_ptr = conf_model[chain->index[ic]];
                //struct ForegroundModel *conf_trial_ptr = conf_trial[chain->index[ic]];
                struct SGWBModel       *sgwb_model_ptr = sgwb_model[chain->index[ic]];
                struct SGWBModel       *sgwb_trial_ptr = sgwb_trial[chain->index[ic]];

                for(int mc=0; mc<10; mc++)
                {
                    noise_instrument_model_mcmc(orbit, data, inst_model_ptr, inst_trial_ptr, conf_model_ptr, sgwb_model_ptr, psd_ptr, chain, flags, ic);
                    //if(flags->confNoise) noise_foreground_model_mcmc(data, inst_model_ptr, conf_model_ptr, conf_trial_ptr, sgwb_model_ptr, psd_ptr, chain, flags, ic);
                    if(flags->sgwbTemplate>=0) noise_sgwb_model_mcmc(data, inst_model_ptr, conf_model_ptr, sgwb_model_ptr, sgwb_trial_ptr, psd_ptr, chain, flags, ic);
                }
            }// end (parallel) loop over chains
            
            //Next section is single threaded. Every thread must get here before continuing
            
            #pragma omp barrier
            
            if(threadID==0)
            {
                noise_ptmcmc(inst_model, chain, flags);
                
                if(step%(flags->NMCMC/10)==0)printf("noise_mcmc at step %i\n",step);
                
                // print chain files
                fprintf(noiseChainFile,"%i %.12g ",step,inst_model[chain->index[0]]->logL);
                print_instrument_state(inst_model[chain->index[0]], noiseChainFile);
                fprintf(noiseChainFile,"\n");

                if(flags->confNoise) 
                {
                    fprintf(foregroundChainFile,"%i %.12g ",step,conf_model[chain->index[0]]->logL);
                    print_foreground_state(conf_model[chain->index[0]], foregroundChainFile);
                    fprintf(foregroundChainFile,"\n");
                }

                if(flags->sgwbTemplate>=0) 
                {
                    fprintf(sgwbChainFile,"%i %.12g ",step,sgwb_model[chain->index[0]]->logL);
                    print_sgwb_state(sgwb_model[chain->index[0]], sgwbChainFile);
                    fprintf(sgwbChainFile,"\n");
                }

                if(step%(flags->NMCMC/10)==0)
                {
                    generate_instrument_noise_model(orbit,inst_model[chain->index[0]]);
                    sprintf(filename,"%s/current_instrument_noise_model.dat",data->dataDir);
                    print_noise_model(inst_model[chain->index[0]]->psd, filename);

                    if(flags->confNoise)
                    {
                        generate_galactic_foreground_model(conf_model[chain->index[0]]);
                        sprintf(filename,"%s/current_foreground_noise_model.dat",data->dataDir);
                        print_noise_model(conf_model[chain->index[0]]->psd, filename);
                    }
                    if(flags->sgwbTemplate>=0) 
                    {
                        generate_sgwb_model(sgwb_model[chain->index[0]]);
                        sprintf(filename,"%s/current_sgwb_noise_model.dat",data->dataDir);
                        print_noise_model(sgwb_model[chain->index[0]]->psd, filename);
                    }
                }
                
                if(step%data->downsample==0 && step/data->downsample < data->Nwave)
                {
                    generate_instrument_noise_model(orbit,inst_model[chain->index[0]]);
                    if(flags->confNoise)
                    {
                        generate_galactic_foreground_model(conf_model[chain->index[0]]);
                        generate_full_covariance_matrix(inst_model[chain->index[0]]->psd,conf_model[chain->index[0]]->psd, data->Nchannel);
                    }
                    if(flags->sgwbTemplate>=0) 
                    {
                        generate_sgwb_model(sgwb_model[chain->index[0]]);
                        generate_full_covariance_matrix(inst_model[chain->index[0]]->psd,sgwb_model[chain->index[0]]->psd, data->Nchannel);
                    }
                    
                    for(int n=0; n<data->NFFT; n++)
                        for(int i=0; i<data->Nchannel; i++)
                            data->S_pow[n][i][step/data->downsample] = inst_model[chain->index[0]]->psd->C[i][i][n];
                }

                step++;
                
                
            }
            //Can't continue MCMC until single thread is finished
            #pragma omp barrier
            
        }// end of MCMC loop
        
    }// End of parallelization
    
    fclose(noiseChainFile);
    if(flags->confNoise)fclose(foregroundChainFile);
    if(flags->sgwbTemplate>=0)fclose(sgwbChainFile);

    generate_instrument_noise_model(orbit,inst_model[chain->index[0]]);
    sprintf(filename,"%s/final_instrument_noise_model.dat",data->dataDir);
    print_noise_model(inst_model[chain->index[0]]->psd, filename);

    if(flags->confNoise)
    {
        generate_galactic_foreground_model(conf_model[chain->index[0]]);
        sprintf(filename,"%s/final_foreground_noise_model.dat",data->dataDir);
        print_noise_model(conf_model[chain->index[0]]->psd, filename);
        generate_full_covariance_matrix(inst_model[chain->index[0]]->psd, conf_model[chain->index[0]]->psd, data->Nchannel);
    }
    if(flags->sgwbTemplate>=0)
    {
        generate_sgwb_model(sgwb_model[chain->index[0]]);
        sprintf(filename,"%s/final_sgwb_noise_model.dat",data->dataDir);
        print_noise_model(sgwb_model[chain->index[0]]->psd, filename);
        generate_full_covariance_matrix(inst_model[chain->index[0]]->psd, sgwb_model[chain->index[0]]->psd, data->Nchannel);
    }
    if(flags->confNoise || flags->sgwbTemplate>=0)
    {
        sprintf(filename,"%s/final_full_noise_model.dat",data->dataDir);
        print_noise_model(inst_model[chain->index[0]]->psd, filename);
    }
    
    print_noise_reconstruction(data, flags);

    sprintf(filename,"%s/whitened_data.dat",data->dataDir);
    // this was here, but it was already added??
    //if(flags->confNoise)generate_full_covariance_matrix(inst_model[chain->index[0]]->psd,conf_model[chain->index[0]]->psd, data->Nchannel);
    print_whitened_data(data, inst_model[chain->index[0]]->psd, filename);

    
    //print total run time
    stop = time(NULL);
    
    printf(" ELAPSED TIME = %g seconds\n",(double)(stop-start));
    
    return 0;
}
