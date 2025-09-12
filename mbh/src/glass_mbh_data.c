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

/**
 @file glass_mbh_data.c
 \brief Data handling functions for MBH module
 */

/*  REQUIRED LIBRARIES  */

#include <glass_utils.h>
#include <glass_noise.h>
#include "glass_mbh.h"

void MBHInjectSimulatedSource(struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Source **inj_vec)
{
    int k;
    
    FILE *fptr;
    FILE *injectionFile;
    FILE *paramFile;
    char filename[1024];
    
    struct TDI *tdi = data->tdi;
    struct Source *inj = malloc(sizeof(struct Source));
    alloc_source(inj, data->N, MBH_MODEL_NP, data->Nchannel);
    
    //reset wavelet basis max and min ranges
    wavelet_pixel_to_index(data->wdm,0,data->lmin,&data->wdm->kmin);
    wavelet_pixel_to_index(data->wdm,0,data->lmax,&data->wdm->kmax);
    
    //recompute fmin and fmax so they align with a bin
    data->fmin = data->lmin*WAVELET_BANDWIDTH;
    data->fmax = data->lmax*WAVELET_BANDWIDTH;
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin+data->NFFT;
    
    fprintf(stdout,"  Frequency layers [%i,%i]\n",data->lmin,data->lmax);
    fprintf(stdout,"  fmin=%lg, fmax=%lg\n",data->fmin,data->fmax);


    int n_inj = 0; //keep track of number of injections
    for(int m=0; m<flags->NINJ; m++)
    {
        injectionFile = fopen(flags->injFile[m],"r");
        if(!injectionFile)
            fprintf(stderr,"Missing injection file %s\n",flags->injFile[m]);
        else
            fprintf(stdout,"Injecting simulated source %s  (%i/%i)\n",flags->injFile[m],m+1, flags->NINJ);

        //count sources in the file
        int N=0;
        while(!feof(injectionFile))
        {
            scan_mbh_source_params(data, inj, injectionFile);
            N++;
        }
        rewind(injectionFile);
        N--;

        for(int n=0; n<N; n++)
        {
            scan_mbh_source_params(data, inj, injectionFile);
            
            //save parameters to file
            sprintf(filename,"%s/injection_parameters_%i.dat",flags->runDir,m);
            if(n==0)paramFile=fopen(filename,"w");
            else    paramFile=fopen(filename,"a");
            print_mbh_source_params(data, inj, paramFile);
            fprintf(paramFile,"\n");
            fclose(paramFile);

            //set analysis trigger time as merger time of injection
            flags->triggerTime = inj->tref;
            
            //generate injection waveform
            mbh_fd_waveform(orbit, data->wdm, data->T, data->t0, inj->params, inj->list, &inj->Nlist, inj->tdi->X, inj->tdi->Y, inj->tdi->Z);

            //inject waveform into data
            for(int j=0; j<inj->Nlist; j++)
            {
                k = inj->list[j];
                if(k>=0 && k<data->N)
                {
                    tdi->X[k] += inj->tdi->X[k];
                    tdi->Y[k] += inj->tdi->Y[k];
                    tdi->Z[k] += inj->tdi->Z[k];
                }
            }
            
            //save injection waveform to file
            sprintf(filename,"%s/data/waveform_injection_%i.dat",flags->runDir,n_inj);
            fptr=fopen(filename,"w");
            for(int j=data->lmin; j<data->lmax; j++)
            {
                double f = j*data->wdm->df;
                for(int i=0; i<data->wdm->NT; i++)
                {
                    double t = i*data->wdm->dt;
                    wavelet_pixel_to_index(data->wdm,i,j,&k);
                    k-=data->wdm->kmin;
                    fprintf(fptr,"%.14e %.14e %.14e %.14e %.14e\n", t, f, inj->tdi->X[k], inj->tdi->Y[k], inj->tdi->Z[k]);
                }
                fprintf(fptr,"\n");
            }
            fclose(fptr);
            
            //save fourier domain injection waveform
            sprintf(filename,"%s/waveform_injection_fourier_%i.dat",data->dataDir,n_inj);
            print_wavelet_fourier_spectra(data, inj->tdi, filename);
            
            //save wavelet domain scaleogram of injection waveform
            sprintf(filename,"%s/data/power_injection_%i.dat",flags->runDir,n_inj);
            fptr=fopen(filename,"w");
            for(int j=data->lmin; j<data->lmax; j++)
            {
                double f = j*data->wdm->df;
                for(int i=0; i<data->wdm->NT; i++)
                {
                    double t = i*data->wdm->dt;
                    wavelet_pixel_to_index(data->wdm,i,j,&k);
                    k-=data->wdm->kmin;
                    fprintf(fptr,"%.14e %.14e %.14e %.14e %.14e\n", t, f, inj->tdi->X[k]*inj->tdi->X[k], inj->tdi->Y[k]*inj->tdi->Y[k], inj->tdi->Z[k]*inj->tdi->Z[k]);
                }
                fprintf(fptr,"\n");
            }
            fclose(fptr);
            
            //Get noise spectrum for data segment
            GetDynamicNoiseModel(data, orbit, flags);

            //Get injected SNR
            if(!flags->quiet)
                fprintf(stdout,"   ...injected SNR=%g\n",snr_wavelet(inj,data->noise));

            
            //Compute fisher information matrix of injection
            if(!flags->quiet)fprintf(stdout,"   ...computing Fisher Information Matrix of injection\n");
            mbh_fisher(orbit, data, inj, data->noise);
            
            if(!flags->quiet)
            {
                fprintf(stdout,"\n Fisher Matrix:\n");
                for(int i=0; i<MBH_MODEL_NP; i++)
                {
                    fprintf(stdout," ");
                    for(int j=0; j<MBH_MODEL_NP; j++)
                    {
                        fprintf(stdout,"%+.2e ", inj->fisher_matrix[i][j]);
                    }
                    fprintf(stdout,"\n");
                }
                
                
                printf("\n Fisher std. errors:\n");
                invert_matrix(inj->fisher_matrix, MBH_MODEL_NP);
                for(int j=0; j<MBH_MODEL_NP; j++)  fprintf(stdout," %.2e\n", sqrt(inj->fisher_matrix[j][j]));
                //exit(1);
            }
            
            if(n_inj<flags->DMAX)
            {
                alloc_source(inj_vec[n_inj], data->N, MBH_MODEL_NP, data->Nchannel);
                copy_source(inj,inj_vec[n_inj]);
            }
            else
            {
                printf("WARNING: number of injections (%i) exceeds size of model (%i)\n",n_inj,flags->DMAX);
                copy_source(inj,inj_vec[flags->DMAX-1]);
            }
            n_inj++;

        }// end loop over source file
        
        fclose(injectionFile);
    }
    
    sprintf(filename,"%s/data/waveform_injection_fourier.dat",flags->runDir);
    print_wavelet_fourier_spectra(data, inj_vec[0]->tdi, filename);
    
    if(!flags->quiet)fprintf(stdout,"================================================\n\n");
}
