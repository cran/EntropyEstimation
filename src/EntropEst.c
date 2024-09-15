//
//  Created by Lijuan Cao and Michael Grabchak on 11/8/14.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include<R.h>



void KlPlugin(int *c1, int *c2, int *h, double *result)
{
    double kp = 0;
    int len = *h;
    double *q1, *q2;
    q1= (double *)R_Calloc(sizeof(double) * len, double);
    q2= (double *)R_Calloc(sizeof(double) * len, double);
  
    int samples1 = 0, samples2 = 0;
    
    for(int i = 0; i < len; i++)
    {
        samples1 +=c1[i];
        samples2 += c2[i];
    }
    
    
    for(int i = 0; i < len; i++)
    {
        q1[i] = c1[i]/((double)samples1);
        q2[i] = c2[i]/((double)samples2);
        if(q2[i]==0)
            q2[i]=1/((double)samples2);
    }
    
    
    for(int i = 0; i< len; i++)
    {
        if(c1[i] ==0)
            continue;
        kp += q1[i] * (log(q1[i])-log(q2[i]));
        
    }
    
    *result = kp;
    free(q1);
    free(q2);
    
}

void KlSharp(int *c1, int *c2, int *h, double *result)
{
    double ks =0;
    int samples1 = 0, samples2 = 0;
    int len = *h;
    
    for(int i = 0; i < len; i++)
    {
        samples1 +=c1[i];
        samples2 += c2[i];
    }
    
    
    for(int i = 0; i< len; i++)
    {
        if(c1[i] ==0)
            continue;
        else{
            double temp = 1;
            double temp2 = 0;
            for(int k = 1; k<=(samples2-c2[i]); k++)
            {
                temp *= 1-((double)c2[i])/(samples2-k+1);
                temp2 += temp/k;
                
            }
            temp = 1;
            for(int k = 1; k<=(samples1-c1[i]); k++)
            {
                temp *= 1- (  ((double)c1[i]) - 1  )/(samples1-k);
                temp2 -= temp/k;
                
            }
            
            ks += temp2 * ((double)c1[i])/samples1;
        }//else
    }
    
    
    *result = ks;
}


void KlSd(int *c1orig, int *c2orig, int *h, double *result)
{
    int len = *h;
    double *g, *c1, *c2;
    g = (double *)R_Calloc(sizeof(double) * (2*(len - 1)), double);
    c1= (double *)R_Calloc(sizeof(double) * len, double);
    c2= (double *)R_Calloc(sizeof(double) * len, double);
    
    int samples1 = 0, samples2 = 0;
    
    for(int i = 0; i< len; i++)
    {
        samples1 +=c1orig[i];
        samples2 += c2orig[i];
        
        c2[i] = c2orig[i];
        c1[i] = c1orig[i];
        if(c1[i]==0 && c2[i] != 0)
            c1[i]=1;
        if(c1[i]!=0 &&c2[i] ==0)
            c2[i] = 1;
    }
    
    int index = 0;
    for(int i = len-1; i>=0; i--)
    {
        if(c1[i] !=0)
        {
            index = i;
            break;
        }
        
    }
    
    for(int i = 0; i < len-1; i++)
    {
        if(c1[i] ==0)
        {
            g[i] = 0;
            g[len-1+i] = 0;
        }
        else
        {
            g[i] = log( c1[i]/((double)c2[i]) ) - log( c1[index]/((double)c2[index]) );
            g[len-1 + i] = -c1[i]*samples2/(((double)c2[i]) * samples1) + c1[index]*samples2/( ( (double)c2[index]) * samples1);
        }
    }
    
    double **Sigma1, **Sigma2;
    Sigma1 = (double**)R_Calloc(sizeof(double*)*(len-1), double);
    for(int i = 0; i < len-1; i++)
        Sigma1[i] = (double *)R_Calloc(sizeof(double)*(len-1), double);

    Sigma2 = (double**)R_Calloc(sizeof(double*)*(len-1), double);
    for(int i = 0; i < len-1; i++)
        Sigma2[i] = (double *)R_Calloc(sizeof(double)*(len-1), double);
    
  
    for(int i = 0; i < len -1; i++)
    {
        for(int j = 0; j < len-1; j++)
        {
            if(i==j)
            {
                Sigma1[i][j] = c1[i]/((double)samples1)*(1-c1[i]/((double)samples1));
                Sigma2[i][j] = c2[i]/((double)samples2)*(1-c2[i]/((double)samples2));
            }
            else
            {
                Sigma1[i][j] = -c1[i]*c1[j]/((double)samples1*samples1);
                Sigma2[i][j] = -c2[i]*c2[j]/((double)samples2*samples2);
            }
        }
    }
    
    double var = 0;
    for(int i = 0; i< len-1; i++)
    {
        for(int j = 0; j< len-1; j++)
        {
            var += g[i]*Sigma1[i][j]*g[j] + g[len-1+i]*Sigma2[i][j]*g[len-1+j];
        }
    }
    

    *result = sqrt(var);
    
    for(int i = 0; i < len-1; i++)
        free(Sigma1[i]);
    
    free(Sigma1);
    
    for(int i = 0; i < len-1; i++)
        free(Sigma2[i]);
    
    free(Sigma2);
    free(g);
    free(c1);
    free(c2);
   
}

void SymSd(int *c1orig, int *c2orig, int *h, double *result)
{
    int len = *h;
    double *g, *c1, *c2;
    g = (double *)R_Calloc(sizeof(double) * (2*(len - 1)), double);
    c1= (double *)R_Calloc(sizeof(double) * len, double);
    c2= (double *)R_Calloc(sizeof(double) * len, double);
    
    int samples1 = 0, samples2 = 0;
    
    for(int i = 0; i< len; i++)
    {
        samples1 +=c1orig[i];
        samples2 += c2orig[i];
        
        c2[i] = c2orig[i];
        c1[i] = c1orig[i];
        if(c1[i]==0 && c2[i] != 0)
            c1[i]=1;
        if(c1[i]!=0 &&c2[i] ==0)
            c2[i] = 1;
    }
    int index = 0;
    for(int i = len-1; i>=0; i--)
    {
        if(c1[i] !=0)
        {
            index = i;
            break;
        }
        
    }
    for(int i = 0; i < len-1; i++)
    {
        if(c1[i] ==0)
        {
            g[i] = 0;
            g[len-1+i] = 0;
        }
        else
        {
            g[i] = .5*( log(c1[i]/((double)c2[i])) - log(c1[index]/((double)c2[index])) ) - .5*(c2[i]*samples1/(((double)c1[i]) * samples2) - c2[index]*samples1/( ((double)c1[index]) * samples2 ) );
            g[len-1 + i] = .5*(log(c2[i]/((double)c1[i])) - log(c2[index]/((double)c1[index])) ) - .5*( c1[i]*samples2/(((double)c2[i])*samples1) - c1[index]*samples2/(((double)c2[index])* samples1)  );
        }
    }
    
    double **Sigma1, **Sigma2;
    Sigma1 = (double **)R_Calloc(sizeof(double*)*(len-1), double);
    for(int i = 0; i < len-1; i++)
        Sigma1[i] = (double *)R_Calloc(sizeof(double)*(len-1), double);
    
    Sigma2 = (double **)R_Calloc(sizeof(double*)*(len-1), double);
    for(int i = 0; i < len-1; i++)
        Sigma2[i] = (double *)R_Calloc(sizeof(double)*(len-1), double);

    for(int i = 0; i < len -1; i++)
    {
        for(int j = 0; j < len-1; j++)
        {
            if(i==j)
            {
                Sigma1[i][j] = c1[i]/((double)samples1)*(1-c1[i]/((double)samples1));
                Sigma2[i][j] = c2[i]/((double)samples2)*(1-c2[i]/((double)samples2));
            }
            else
            {
                Sigma1[i][j] = -c1[i]*c1[j]/((double)samples1*samples1);
                Sigma2[i][j] = -c2[i]*c2[j]/((double)samples2*samples2);
            }
        }
    }
    
    double var = 0;
    for(int i = 0; i < len-1; i++)
    {
        for(int j = 0; j < len-1; j++)
        {
            var += g[i]*Sigma1[i][j]*g[j] + g[len-1+i]*Sigma2[i][j]*g[len-1+j];
        }
    }
    

    *result = sqrt(var);
    
    for(int i = 0; i < len-1; i++)
        free(Sigma1[i]);
    
    free(Sigma1);
    
    for(int i = 0; i < len-1; i++)
        free(Sigma2[i]);
    
    free(Sigma2);
    free(g);
    free(c1);
    free(c2);

   
}

void EntropySharp(int *c1, int *h, double *result)
{
    double ent =0;
    int samples1 = 0;
    int len = *h;
    
    for(int i = 0; i < len; i++)
    {
        samples1 +=c1[i];
    }
    
    
    for(int i = 0; i< len; i++)
    {
        if(c1[i] ==0)
            continue;
        else
        {
            double temp = 1;
            double temp2 = 0;
            for(int k = 1; k<=(samples1-c1[i]); k++)
            {
                temp *= 1- (  ((double)c1[i]) - 1  )/(samples1-k);
                temp2 += temp/k;
                
            }
            
            ent += temp2 * ((double)c1[i])/samples1;
        }//else
    }
    
    
    *result = ent;
}


void RenyiEqEntropySharp(int *c1, int *h, double *r, double *result)
{
    double ent =0;
    int samples1 = 0;
    int len = *h;
    
    for(int i = 0; i < len; i++)
    {
        samples1 +=c1[i];
    }
    
    double w;
    
    for(int i = 0; i< len; i++)
    {
        if(c1[i] ==0)
            continue;
        else
        {
            w=1;
            double temp = 1;
            double temp2 = 0;
            for(int k = 1; k<=(samples1-c1[i]); k++)
            {
                w *= (1- *r/k);
                temp *= 1- (  ((double)c1[i]) - 1  )/(samples1-k);
                temp2 += w * temp;
                
            }
            
            ent += temp2 * ((double)c1[i])/samples1;
        }//else
    }
    
    
    *result = 1+ent;
}


void EntropySd(int *c1orig, int *h, double *result)
{
    
    int len = *h;
    double *g, *c1;
    g = (double *)R_Calloc(sizeof(double) * (len - 1), double);
    c1= (double *)R_Calloc(sizeof(double) * len, double);
  
    
    int samples1 = 0;
    
    for(int i = 0; i< len; i++)
    {
        samples1 +=c1orig[i];
        c1[i] = c1orig[i];
    }
    
    int index = 0;
    for(int i = len-1; i>=0; i--)
    {
        if(c1[i] !=0)
        {
            index = i;
            break;
        }
        
    }
    
    for(int i = 0; i < len-1; i++)
    {
        if(c1[i] ==0)
        {
            g[i] = 0;
        }
        else
        {
            g[i] = log( c1[i]/((double)c1[index]) );
        }
    }
    double **Sigma1;
    Sigma1 = (double **)R_Calloc(sizeof(double*)*(len-1), double);
    for(int i = 0; i < len-1; i++)
        Sigma1[i] = (double *)R_Calloc(sizeof(double)*(len-1), double);

    for(int i = 0; i < len -1; i++)
    {
        for(int j = 0; j < len-1; j++)
        {
            if(i==j)
            {
                Sigma1[i][j] = c1[i]/((double)samples1)*(1-c1[i]/((double)samples1));
            }
            else
            {
                Sigma1[i][j] = -c1[i]*c1[j]/((double)samples1*samples1);
            }
        }
    }
    
    double var = 0;
    for(int i = 0; i< len-1; i++)
    {
        for(int j = 0; j< len-1; j++)
        {
            var += g[i]*Sigma1[i][j]*g[j];
        }
    }
    
    
    *result = sqrt(var);
    
    for(int i = 0; i < len-1; i++)
        free(Sigma1[i]);
    
    free(Sigma1);
    free(g);
    free(c1);
}

void RenyiEqSd(int *c1orig, int *h, double *r, double *result)
{
    int len = *h;
    double *g, *c1;
    g = (double *)R_Calloc(sizeof(double) * (len - 1), double);
    c1= (double *)R_Calloc(sizeof(double) * len, double);

    
    int samples1 = 0;
    
    for(int i = 0; i< len; i++)
    {
        samples1 +=c1orig[i];
        c1[i] = c1orig[i];
    }
    
    int index = 0;
    for(int i = len-1; i>=0; i--)
    {
        if(c1[i] !=0)
        {
            index = i;
            break;
        }
        
    }
    
    for(int i = 0; i < len-1; i++)
    {
        if(c1[i] ==0)
        {
            g[i] = 0;
        }
        else
        {
            g[i] = *r *( pow(c1[i]/((double)samples1), *r-1) + pow(c1[index]/((double)samples1), *r-1) );
        }
    }
    
    double **Sigma1;
    Sigma1 = (double **)R_Calloc(sizeof(double*)*(len-1), double);
    for(int i = 0; i < len-1; i++)
        Sigma1[i] = (double *)R_Calloc(sizeof(double)*(len-1), double);
    
    for(int i = 0; i < len -1; i++)
    {
        for(int j = 0; j < len-1; j++)
        {
            if(i==j)
            {
                Sigma1[i][j] = c1[i]/((double)samples1)*(1-c1[i]/((double)samples1));
            }
            else
            {
                Sigma1[i][j] = -c1[i]*c1[j]/((double)samples1*samples1);
            }
        }
    }
    
    double var = 0;
    for(int i = 0; i< len-1; i++)
    {
        for(int j = 0; j< len-1; j++)
        {
            var += g[i]*Sigma1[i][j]*g[j];
        }
    }
    
    
    *result = sqrt(var);
    
    for(int i = 0; i < len-1; i++)
        free(Sigma1[i]);
    
    free(Sigma1);
    free(g);
    free(c1);
    
}

void GenSimpSharp(int *c1, int *h, int *r, double *result)
{
    double ent =0;
    int samples1 = 0;
    int len = *h;
    
    for(int i = 0; i < len; i++)
    {
        samples1 +=c1[i];
    }
    
    
    for(int i = 0; i< len; i++)
    {
        if(c1[i] ==0)
            continue;
        else
        {
            double temp = 1;
            for(int k = 1; k<= *r; k++)
            {
                temp *= 1- (  ((double)c1[i]) - 1  )/(samples1-k);
            }
            
            ent += temp * ((double)c1[i])/samples1;
        }//else
    }
    
    
    *result = ent;
}

void MISd(int *c1orig, int *h, double *g, double *result)
{
    int len = *h;
    double *c1;
    c1= (double *)R_Calloc(sizeof(double) * len, double);
    
    int samples1 = 0;
    
    for(int i = 0; i< len; i++)
    {
        samples1 +=c1orig[i];
        c1[i] = c1orig[i];
    }
    
    double **Sigma1;
    Sigma1 = (double **)R_Calloc(sizeof(double*)*(len-1), double);
    for(int i = 0; i < len-1; i++)
        Sigma1[i] = (double *)R_Calloc(sizeof(double)*(len-1), double);
    
    for(int i = 0; i < len -1; i++)
    {
        for(int j = 0; j < len-1; j++)
        {
            if(i==j)
            {
                Sigma1[i][j] = c1[i]/((double)samples1)*(1-c1[i]/((double)samples1));
            }
            else
            {
                Sigma1[i][j] = -c1[i]*c1[j]/((double)samples1*samples1);
            }
        }
    }
    
    double var = 0;
    for(int i = 0; i< len-1; i++)
    {
        for(int j = 0; j< len-1; j++)
        {
            var += g[i]*Sigma1[i][j]*g[j];
        }
    }
    
    
    *result = sqrt(var);
    
    for(int i = 0; i < len-1; i++)
        free(Sigma1[i]);
    
    free(Sigma1);
    free(c1);
    
}

void GenSimpSd(int *c1orig, int *h, int *r, double *result)
{
    int len = *h;
    double *g, *c1;
    g = (double *)R_Calloc(sizeof(double) * (len - 1), double);
    c1= (double *)R_Calloc(sizeof(double) * len, double);
    
    int samples1 = 0;
    
    for(int i = 0; i< len; i++)
    {
        samples1 +=c1orig[i];
        c1[i] = c1orig[i];
    }
    
    int index = 0;
    for(int i = len-1; i>=0; i--)
    {
        if(c1[i] !=0)
        {
            index = i;
            break;
        }
        
    }
    
    for(int i = 0; i < len-1; i++)
    {
        if(c1[i] ==0)
        {
            g[i] = 0;
        }
        else
        {
            g[i] = pow((1 - c1[i]/samples1), *r) - *r * (c1[i]/samples1)*pow((1 - c1[i]/samples1), *r - 1) - pow((1 - c1[index]/samples1), *r) + *r *(c1[index]/samples1)* pow((1 - c1[index]/samples1), *r - 1);
        }
    }
    
    double **Sigma1;
    Sigma1 = (double **)R_Calloc(sizeof(double*)*(len-1), double);
    for(int i = 0; i < len-1; i++)
        Sigma1[i] = (double *)R_Calloc(sizeof(double)*(len-1), double);
    
    for(int i = 0; i < len -1; i++)
    {
        for(int j = 0; j < len-1; j++)
        {
            if(i==j)
            {
                Sigma1[i][j] = c1[i]/((double)samples1)*(1-c1[i]/((double)samples1));
            }
            else
            {
                Sigma1[i][j] = -c1[i]*c1[j]/((double)samples1*samples1);
            }
        }
    }
    
    double var = 0;
    for(int i = 0; i< len-1; i++)
    {
        for(int j = 0; j< len-1; j++)
        {
            var += g[i]*Sigma1[i][j]*g[j];
        }
    }
    
    
    *result = sqrt(var);
    
    
    for(int i = 0; i < len-1; i++)
        free(Sigma1[i]);
    
    free(Sigma1);
    free(g);
    free(c1);

    
}


