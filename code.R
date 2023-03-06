sorgente <- function(stoc,n_it,Burnin,Thin) { 
  ##Function to be used together with the dataset provided
  ##Input:
  #stoc=1 to apply stochastic method, stoc=0 to apply deterministic method
  #n_it=Total number of iterations of the MCMC (for stoc. method)
  #Burnin=How much burnin to apply (for stoc. method)
  #Thin=How much thin to apply (for stoc. method)
  if(stoc==0){ #Det. method
    library(stats4)
    I_o_t <- function(k,epsilon) { #epsilon stands for parameter I0
      Iot=numeric(n)
      Iot[1]=epsilon
      for (i in 2:n) {
        Iot[i]=Iot[i-1]+(k*Iot[i-1]/(I_tot_t[i-1]+(k-1)*Iot[i-1]))*Y_tot_t[i]-
          (Iot[i-1]/I_tot_t[i-1])*diff_R[i]
      }
      return(Iot)
    }
    l <- function(k,epsilon) { 
      Io_t=I_o_t(k,epsilon)
      Y0=numeric(n)
      p=numeric(n)
      loglik=0
      for (t in 2:n) {
        p[t]=k*Io_t[t-1]/(I_tot_t[t-1]+(k-1)*Io_t[t-1])
        Y0[t]=p[t]*Y_tot_t[t]
        loglik=loglik+lchoose(Y0[t],z[t])+lchoose(Y_tot_t[t]-Y0[t],n_t[t]-z[t])
      }
      return(-loglik)
    }
    est<-mle(minuslogl=l,method = "L-BFGS-B",lower=c(0,1),
             start = list(k = 2, epsilon = 10))
    sol=matrix(0,2,5) 
    #the matrix sol will contain:
    #in the first row -> estimate of k, 95% asymptotic C.I., 95% C.I. from parametric bootstrap
    #in the second row -> same for parameter I0 (called also epsilon)
    sol[,1]=est@coef
    conf=confint(est)
    sol[1,2:3]=conf[1,]
    sol[2,2:3]=conf[2,]
    
    ## parametric bootstrap
    k=est@coef[1]
    eps=est@coef[2]
    Io_t=I_o_t(k,epsilon=eps)
    inc=numeric(n)
    Yom=numeric(n)
    for (t in 2:n) {
      inc[t]=k*Io_t[t-1]/(I_tot_t[t-1]+(k-1)*Io_t[t-1])
      Yom[t]=inc[t]*Y_tot_t[t]
    }
    
    x_axis=c(2,16,30,44,58)
    N=1000
    zz=matrix(0,ncol=n,nrow=N)
    for(day in x_axis){
      zz[,day]<-rhyper(N,Yom[day],Y_tot_t[day]-Yom[day],n_t[day])
    }
    
    kk=numeric(N)
    epss=numeric(N)
    library(extraDistr)
    for(i in 1:N){
      z_use=zz[i,]
      l_bis <- function(k,epsilon) { 
        Io_t=I_o_t(k,epsilon)
        Y0=numeric(n)
        p=numeric(n)
        loglik=0
        for (t in 2:n) {
          p[t]=k*Io_t[t-1]/(I_tot_t[t-1]+(k-1)*Io_t[t-1])
          Y0[t]=p[t]*Y_tot_t[t]
          loglik=loglik+lchoose(Y0[t],z_use[t])+lchoose(Y_tot_t[t]-Y0[t],n_t[t]-z_use[t])
        }
        return(-loglik)
      }
      est<-mle(minuslogl=l_bis,method = "L-BFGS-B",lower=c(0,1),
               start = list(k = 1, epsilon = 3))
      kk[i]=est@coef[1]
      epss[i]=est@coef[2]
    }
    #C.I. 95% from parametric bootstrap:
    sol[1,4]=quantile(kk,0.025)
    sol[1,5]=quantile(kk,0.975)
    sol[2,4]=quantile(epss,0.025)
    sol[2,5]=quantile(epss,0.975)
    
  }else{
    if(stoc==1){ #Stoc. method
      iterazioni=n_it
      burnin=Burnin
      thin=Thin
      N = floor((iterazioni-burnin)/thin)
      #the following variables were used to tune the acceptance rates of the parameters
      eps_acc=0 #epsilon stands for I0 (memorised in position 1)
      k_acc=0
      y_acc=numeric(n)
      i_acc=numeric(n)
      # To memorize the samples
      k = matrix(0,ncol=1,nrow=N) 
      Y = matrix(0,ncol=n,nrow=N)
      I = matrix(0,ncol=n,nrow=N)
      #starting values, and setting them as current values of the MCMC
      k_start=2
      eps_start=10
      y_corr=z
      rho_corr=numeric(n)
      delta_i_corr=y_corr-rho_corr
      delta_i_corr[1]=eps_start
      i_corr=cumsum(delta_i_corr)
      k_corr=k_start
      #For the uniform proposal distributions:
      Nn=c(rep(20,5),rep(30,4),rep(50,7),rep(70,3),rep(80,10),rep(100,10),rep(80,5),rep(30,n-44))
      Nrho=c(rep(20,5),rep(20,4),rep(20,7),rep(50,3),rep(70,10),rep(130,11),rep(150,n-40))
      N_eps=20 
      
      ##### Let us define some functions:
      #full-cond. distribution of y_t (t from 2 to n-1)
      fc_y_total <- function(t,y) {  #t=time, y=y_corr (vector)
        yt=y[t]
        y_d=Y_tot_t[t]-yt
        rt=i_corr[t-1]+yt-i_corr[t]
        somma=lchoose(yt,z[t])+lchoose(y_d,n_t[t]-z[t])+lchoose(Y_tot_t[t],yt)+
          yt*log(k_corr*i_corr[t-1])+y_d*log(I_tot_t[t-1]-i_corr[t-1])-
          Y_tot_t[t]*log(I_tot_t[t-1]+(k_corr-1)*i_corr[t-1])-
          lchoose(Y_tot_t[t],n_t[t])+lchoose(i_corr[t-1],rt)+
          lchoose(I_tot_t[t-1]-i_corr[t-1],diff_R[t]-rt)-lchoose(I_tot_t[t-1],diff_R[t])
        if((exp(somma)<0)|((exp(somma)>1))){
          print('err fc_y_total')
        }
        return(somma)
      }
      
      #full-cond. distribution of y_n (last y_t)
      fc_y_total_final <- function(t,y) {  #t=time, y=y_corr (vector)
        yt=y[t]
        y_d=Y_tot_t[t]-yt
        somma=lchoose(yt,z[t])+lchoose(y_d,n_t[t]-z[t])+lchoose(Y_tot_t[t],yt)+
          yt*log(k_corr*i_corr[t-1])+y_d*log(I_tot_t[t-1]-i_corr[t-1])-
          Y_tot_t[t]*log(I_tot_t[t-1]+(k_corr-1)*i_corr[t-1])-
          lchoose(Y_tot_t[t],n_t[t])
        if((exp(somma)<0)|((exp(somma)>1))){
          print('err fc_y_total_final')
        }
        return(somma)
      }
      
      #full-cond. distribution of i_t (t from 2 to n-2)
      fc_rho_total <- function(t,i) {  #t=time, i=i_corr (vector)
        it=i[t]
        itm1=i[t-1]
        itp1=i[t+1]
        rt=itm1+y_corr[t]-it
        rtp1=it+y_corr[t+1]-itp1
        somma=lchoose(itm1,rt)+lchoose(I_tot_t[t-1]-itm1,diff_R[t]-rt)+
          lchoose(it,rtp1)+
          lchoose(I_tot_t[t]-it,diff_R[t+1]-rtp1)+
          lchoose(Y_tot_t[t+1],y_corr[t+1])+
          y_corr[t+1]*log(k_corr*it)+
          (Y_tot_t[t+1]-y_corr[t+1])*log(I_tot_t[t]-it)-
          Y_tot_t[t+1]*log(I_tot_t[t]+(k_corr-1)*it)-
          lchoose(I_tot_t[t-1],diff_R[t])-lchoose(I_tot_t[t],diff_R[t+1])
        if((exp(somma)<0)|((exp(somma)>1))){
          print('err fc_rho_total')
        }
        return(somma)
      }
      
      #full-cond. distribution of i_{n-1} (last i_t)
      fc_rho_total_final <- function(t,i) {  #t=time, i=i_corr (vector)
        it=i[t]
        itm1=i[t-1]
        rt=itm1+y_corr[t]-it
        somma=lchoose(itm1,rt)+lchoose(I_tot_t[t-1]-itm1,diff_R[t]-rt)+
          lchoose(Y_tot_t[t+1],y_corr[t+1])+
          y_corr[t+1]*log(k_corr*it)+
          (Y_tot_t[t+1]-y_corr[t+1])*log(I_tot_t[t]-it)-
          Y_tot_t[t+1]*log(I_tot_t[t]+(k_corr-1)*it)-lchoose(I_tot_t[t-1],diff_R[t])
        if((exp(somma)<0)|((exp(somma)>1))){
          print('err fc_rho_total_final')
        }
        return(somma)
      }
      
      #full-cond. distribution of i_0 (first i_t, called also epsilon):
      fc_eps <- function(i,y) {  # i e y are vectors
        y1=y[2]
        i0=i[1]
        i1=i[2]
        r1=i0+y1-i1
        somma=lchoose(Y_tot_t[2],y1)+y1*log(k_corr*i0)+
          (Y_tot_t[2]-y1)*log(I_tot_t[1]-i0)-
          Y_tot_t[2]*log(I_tot_t[1]+(k_corr-1)*i0)+
          lchoose(i0,r1)+lchoose(I_tot_t[1]-i0,diff_R[2]-r1)-
          lchoose(I_tot_t[1],diff_R[2])
        if((exp(somma)<0)|((exp(somma)>1))){
          print('err fc_eps')
        }
        return(somma)
      }
      
      #function for the proposal of y_t:
      proposta_y <- function(t,NN,bound_inf,bound_sup,yy) {  
        if((yy[t]>bound_inf+NN) & (yy[t]<bound_sup-NN)){
          #we are far from the two bounds
          y_prop=rdunif(1,yy[t]-NN,yy[t]+NN)
          proposta_num=1/(2*NN+1)
          proposta_den=1/(2*NN+1)
        }else{
          if((yy[t]>bound_inf+NN) & (yy[t]>=bound_sup-NN)){
            #we are far from the lower bound, but close to the upper bound
            y_prop=rdunif(1,yy[t]-NN,bound_sup)
            proposta_num=1/(bound_sup-yy[t]+NN+1)
            proposta_den=1/(bound_sup-y_prop+NN+1)
          }else{
            if((yy[t]<=bound_inf+NN)&(yy[t]<bound_sup-NN)){
              #we are far from the upper bound, but close to the lower bound
              y_prop=rdunif(1,bound_inf,yy[t]+NN)
              proposta_num=1/(yy[t]+NN-bound_inf+1)
              proposta_den=1/(y_prop+NN-bound_inf+1)
            }else{#we are close to both bounds
              y_prop=rdunif(1,bound_inf,bound_sup)
              proposta_num=1/(bound_sup-bound_inf+1)
              proposta_den=1/(bound_sup-bound_inf+1)
            }
          }
        }
        proposta=c(y_prop,proposta_num,proposta_den)
        return(proposta)
      }
      
      #function for the proposal of i_t: 
      proposta_i <- function(t,N_rho,bound_inf,bound_sup,ii) {  
        if((ii[t]>bound_inf+N_rho) & (ii[t]<bound_sup-N_rho)){
          #we are far from the two bounds
          i_prop=rdunif(1,ii[t]-N_rho,ii[t]+N_rho)
          proposta_num=1/(2*N_rho+1)
          proposta_den=1/(2*N_rho+1)
        }else{
          if((ii[t]>bound_inf+N_rho) & (ii[t]>=bound_sup-N_rho)){
            #we are far from the lower bound, but close to the upper bound 
            i_prop=rdunif(1,ii[t]-N_rho,bound_sup)
            proposta_num=1/(bound_sup-ii[t]+N_rho+1)
            proposta_den=1/(bound_sup-i_prop+N_rho+1)
          }else{
            if((ii[t]<=bound_inf+N_rho)&(ii[t]<bound_sup-N_rho)){
              #we are far from the upper bound, but close to the lower bound
              i_prop=rdunif(1,bound_inf,ii[t]+N_rho)
              proposta_num=1/(ii[t]+N_rho-bound_inf+1)
              proposta_den=1/(i_prop+N_rho-bound_inf+1)
            }else{#we are close to both bounds
              i_prop=rdunif(1,bound_inf,bound_sup)
              proposta_num=1/(bound_sup-bound_inf+1)
              proposta_den=1/(bound_sup-bound_inf+1)
            }
          }
        }
        prop_i=c(i_prop,proposta_num,proposta_den)
        return(prop_i)
      }
      
      #function to compute the M.H: alpha rate for parameter k:
      accetto_k_o_no <- function(k_prop,k_corr,y_corr,i_corr) {  
        LogNum=sum(y_corr[2:n]*log(k_prop*i_corr[1:(n-1)])-
                     Y_tot_t[2:n]*log(I_tot_t[1:(n-1)]+(k_prop-1)*i_corr[1:(n-1)]))
        LogDen=sum(y_corr[2:n]*log(k_corr*i_corr[1:(n-1)])-
                     Y_tot_t[2:n]*log(I_tot_t[1:(n-1)]+(k_corr-1)*i_corr[1:(n-1)]))
        LogNum=LogNum+log(k_prop)-0.9*log(k_prop)-0.1*k_prop
        LogDen=LogDen+log(k_corr)-0.9*log(k_corr)-0.1*k_corr
        alpha = min(1,exp(LogNum - LogDen))
        return(alpha)
      }
      
      #function for the proposal of epsilon (I0):
      proposta_eps <- function(eps_corr,b_inf,b_sup) {  
        if((eps_corr>b_inf+N_eps) & (eps_corr<b_sup-N_eps)){
          #we are far from the two bounds
          eps_prop=rdunif(1,eps_corr-N_eps,eps_corr+N_eps)
          proposta_num=1/(2*N_eps+1)
          proposta_den=1/(2*N_eps+1)
        }else{
          if((eps_corr>b_inf+N_eps) & (eps_corr>=b_sup-N_eps)){
            #we are far from the lower bound, but close to the upper bound 
            eps_prop=rdunif(1,eps_corr-N_eps,b_sup)
            proposta_num=1/(b_sup-eps_corr+N_eps+1)
            proposta_den=1/(b_sup-eps_prop+N_eps+1)
          }else{
            if((eps_corr<=b_inf+N_eps) & (eps_corr<b_sup-N_eps)){
              #we are far from the upper bound, but close to the lower bound
              eps_prop=rdunif(1,b_inf,eps_corr+N_eps)
              proposta_num=1/(eps_corr+N_eps-b_inf+1)
              proposta_den=1/(eps_prop+N_eps-b_inf+1)
            }
            else{#we are close to both bounds
              eps_prop=rdunif(1,b_inf,b_sup)
              proposta_num=1/(b_sup-b_inf+1)
              proposta_den=1/(b_sup-b_inf+1)
            }
          }
        }
        prop_eps=c(eps_prop,proposta_num,proposta_den)
        return(prop_eps)
      }
      #####
      
      M=burnin
      for(iter1 in 1:N)
      {
        for(iter2 in 1:M)
        {
          ################################# sample I0 (epsilon)
          b_inf=max(1,i_corr[2]-y_corr[2])
          b_sup=min(I_tot_t[1],i_corr[2]-y_corr[2]+diff_R[2])
          prop_I0=proposta_eps(i_corr[1],b_inf,b_sup)
          eps_prop=prop_I0[1]
          proposta_num=prop_I0[2]
          proposta_den=prop_I0[3]
          i_use=c(eps_prop,i_corr[2:n])
          LogNum=fc_eps(i_use,y_corr)
          LogDen=fc_eps(i_corr,y_corr)
          LogNum=LogNum+log(1/1000)+log(proposta_num)
          LogDen=LogDen+log(1/1000)+log(proposta_den)
          Alpha = min(1,exp(LogNum - LogDen))
          u_eps = runif(1,0,1)
          if(u_eps<Alpha) {
            i_corr[1]=eps_prop
            eps_acc=eps_acc+1 
          }
          #################################
          
          ############################# sample y and i:
          for(t in 2:(n-1)){
            ########################################################## y
            bound_inf=max(z[t],i_corr[t]-i_corr[t-1],
                          diff_R[t]+i_corr[t]-I_tot_t[t-1])
            bound_sup=min(Y_tot_t[t]-n_t[t]+z[t],
                          diff_R[t]+i_corr[t]-i_corr[t-1],i_corr[t])
            prop=proposta_y(t,Nn[t],bound_inf,bound_sup,y_corr)
            y_prop=prop[1]
            proposta_num=prop[2]
            proposta_den=prop[3]
            y_use=c(y_corr[1:(t-1)],y_prop,y_corr[(t+1):n])
            LogNum=fc_y_total(t,y_use)+log(proposta_num)
            LogDen=fc_y_total(t,y_corr)+log(proposta_den)
            Alpha = min(1,exp(LogNum - LogDen))
            u_y = runif(1,0,1)
            if(u_y<Alpha) {
              y_corr[t] = y_prop
              y_acc[t]=y_acc[t]+1
            }
            ###################################################### i
            if(t==(n-1)){ #Last i (time=n-1)
              bound_inf=max(i_corr[t-1]+y_corr[t]-diff_R[t],y_corr[t])
              bound_sup=min(I_tot_t[t],i_corr[t-1]+y_corr[t],
                            I_tot_t[t-1]-diff_R[t]+y_corr[t])
              prop=proposta_i(t,Nrho[t],bound_inf,bound_sup,i_corr)
              i_prop=prop[1]
              proposta_num=prop[2]
              proposta_den=prop[3]
              i_use=c(i_corr[1:(t-1)],i_prop,i_corr[(t+1):n])
              full1=fc_rho_total_final(t,i_use)
              full2=fc_rho_total_final(t,i_corr)
            }else{ # i_t (t from 2 to n-2)
              bound_inf=max(i_corr[t-1]+y_corr[t]-diff_R[t],y_corr[t],
                            i_corr[t+1]-y_corr[t+1])
              bound_sup=min(I_tot_t[t],i_corr[t-1]+y_corr[t],
                            I_tot_t[t-1]-diff_R[t]+y_corr[t],
                            i_corr[t+1]-y_corr[t+1]+diff_R[t+1])
              prop=proposta_i(t,Nrho[t],bound_inf,bound_sup,i_corr)
              i_prop=prop[1]
              proposta_num=prop[2]
              proposta_den=prop[3]
              i_use=c(i_corr[1:(t-1)],i_prop,i_corr[(t+1):n])
              full1=fc_rho_total(t,i_use)
              full2=fc_rho_total(t,i_corr)
            }
            LogNum=full1+log(proposta_num)
            LogDen=full2+log(proposta_den)
            Alpha = min(1,exp(LogNum - LogDen))
            u_i = runif(1,0,1)
            if(u_i<Alpha) {
              i_corr[t] = i_prop
              i_acc[t]=i_acc[t]+1
            }
          }
          ################################ Last y (time=n) 
          t=n
          bound_inf=z[t]
          bound_sup=Y_tot_t[t]-n_t[t]+z[t]
          prop=proposta_y(t,Nn[t],bound_inf,bound_sup,y_corr)
          y_prop=prop[1]
          proposta_num=prop[2]
          proposta_den=prop[3]
          y_use=c(y_corr[1:(t-1)],y_prop)
          LogNum=fc_y_total_final(t,y_use)+log(proposta_num)
          LogDen=fc_y_total_final(t,y_corr)+log(proposta_den)
          Alpha = min(1,exp(LogNum - LogDen))
          u_y = runif(1,0,1)
          if(u_y<Alpha) {
            y_corr[t] = y_prop
            y_acc[t]=y_acc[t]+1
          }
          ########################################################### sample k   
          k_prop=exp(rnorm(1,log(k_corr),0.02))
          Alpha=accetto_k_o_no(k_prop,k_corr,y_corr,i_corr)
          u_k = runif(1,0,1)
          if(u_k<Alpha) {
            k_corr = k_prop
            k_acc=k_acc+1 
          }
        }
        M=thin
        #After the burnin and every 'thin' samples, the samples are memorised
        k[iter1,]=k_corr
        Y[iter1,]=y_corr
        I[iter1,]=i_corr
      }
      sol=cbind(k,Y,I)
      #sol will have a number of rows equal to N
      #The first column contains the samples of k
      #Columns from 2 to n+1 contain the samples of y_t 
      #(we are interested in columns 3:n+1)
      #Columns from n+2 to 2n+1 contain the samples of i_t 
      #(we are interested in columns (n+2):(2n), where column n+2 contains I0)
    }
  }
  return(sol)
}
