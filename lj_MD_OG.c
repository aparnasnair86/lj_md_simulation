// Code for performing Molecular Dynamic Simulation - June 2011

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#define eqmdsteps 50000
#define datamdsteps 100000
#define dt 0.002

long int lround(double x);

int main(int argc, char *argv[]){
	int Num_par,mdstep, atomdone, atomnotdone, N;
	double Rho, lx, ly, lz, xij, yij, zij, rsq, rmindistsq, mindist,rand, rand1, rand2, rand3;
	double *x, *y, *z, *vx, *vy, *vz, *Fx, *Fy, *Fz, *rmin;
	double Fxpar, Fypar, Fzpar;
	double UljFT, UljST,Ulj, Forcexlj, Forceylj, Forcezlj;
	double Uch, Forcexch, Forceych, Forcezch;
	double T, Trt, r, F;

	double xcum, ycum, zcum, vxcum, vycum, vzcum;
	double rc1, rc1sq, rc2, rc2sq;
	// Variables used while initialization
	
	int i, j, k, segnum, num,n;
	double *xdum, ydum, zdum, b0;
	FILE *fpconf, *fpw;


	//variables used while integration
	double *xold, *yold, *zold, *vxold, *vyold, *vzold;
	double *Fxold, *Fyold, *Fzold, pe, ke; 
	const gsl_rng_type *P;
	gsl_rng *m;

	gsl_rng_env_setup();

	P = gsl_rng_default;
	m = gsl_rng_alloc (P);

	//Input Number of particles and density and Temperature //
	if(argc<1+2){
	   printf("Error!!!\n");
	   printf("Format: %s <No.of particles> <Temperature>\n",argv[0]);
   	   exit(1);
	}
	
	num = atoi(argv[1]);
	T = atof(argv[2]);

 	
	xdum = calloc(num, sizeof(double));
	lx = 8.0;
	ly = 8.0;
	lz = 8.0;
	//lx = (double)b0*segnum;

	//fprintf(stderr, "%d\t%f\t%f\n",Num_par, Rho, T);
	x = calloc(num+1, sizeof(double));
	y = calloc(num+1, sizeof(double));
	z = calloc(num+1, sizeof(double));

	xold = calloc(num+1, sizeof(double));
	yold = calloc(num+1, sizeof(double));
	zold = calloc(num+1, sizeof(double));
	
	vx = calloc(num+1, sizeof(double));
	vy = calloc(num+1, sizeof(double));
	vz = calloc(num+1, sizeof(double));

	vxold = calloc(num+1, sizeof(double));
	vyold = calloc(num+1, sizeof(double));
	vzold = calloc(num+1, sizeof(double));
	
	Fx = calloc(num+1, sizeof(double));
	Fy = calloc(num+1, sizeof(double));
	Fz = calloc(num+1, sizeof(double));

	Fxold = calloc(num+1, sizeof(double));
	Fyold = calloc(num+1, sizeof(double));
	Fzold = calloc(num+1, sizeof(double));

	rmin = calloc(num+1, sizeof(double));
	
	//fprintf(stderr, "%f\t%f\t%f\n", lx,ly,lz); 

	/* begin JAD */
	//fpconf = fopen("Init_125_par.dat","r");
	fpconf = fopen("Init_125_par_new.dat","r");
	/* end JAD */

	fpw = fopen("Conf_out.dat","w");

	if(fpconf == NULL){
	  printf("File not opened\n");
	  exit(1);
	}/*else{
	   printf("Opened\n");
	}*/
	rc1 = 1.5; rc1sq = rc1*rc1;
	rc2 = 1.1225; rc2sq = rc2*rc2;
	
//	b0 = 1.0;

//****************Initialize the system with Num_par no.of particles*******************//

	fscanf(fpconf,"%d\n\n", &n);
	//	printf("%d\n", n);
	for(i=1;i<=num;i++){
	   fscanf(fpconf,"%lf\t%lf\t%lf\n", &x[i], &y[i], &z[i]);
	   x[i] = x[i]; //- (0.5*lx);
	   y[i] = y[i]; //- (0.5*lx);
	   z[i] = z[i]; //- (0.5*lx);

	//  printf("Si %lf\t%lf\t%lf\n", x[i], y[i], z[i]);
	}


	//Initialize velocities of the system from Gaussian distribution    //
	Trt = pow(T,1.0/2.0);
	//printf("%d\n", N);

	/* begin JAD */
	vxcum = 0.0;
        vycum = 0.0;
        vzcum = 0.0;
	/* end JAD */

	for(i=1;i<=num;i++){
	   vx[i]=gsl_ran_gaussian(m, Trt);
	   vy[i]=gsl_ran_gaussian(m, Trt);
	   vz[i]=gsl_ran_gaussian(m, Trt);
    	 //  printf("%f\t%f\t%f\n", vx[i],vy[i],vz[i])
	   
	//COM Velocity//
	  
	/* begin JAD */	
	   //vxcum = 0.0;
	   //vycum = 0.0;
	   //vzcum = 0.0;
	/* end JAD */

	   vxcum = vxcum + vx[i];
	   vycum = vycum + vy[i];
	   vzcum = vzcum + vz[i];
	}
	vxcum = vxcum/num;
     	vycum = vycum/num;
	vzcum = vzcum/num;

	for(i=1;i<=num;i++){	   	
  	   vx[i] = vx[i] - vxcum;
	   vy[i] = vy[i] - vycum;
	   vz[i] = vz[i] - vzcum;
	}
// ************************* Initialization done ********************************************* //

	//APARNA----------------------ZERO STEP OF MD-------------------------//

	/*for(i=1;i<=num;i++){
	    for(j=1; j<=num; j++){
	       if(i!=j){
		  xij = (x[i] - x[j]) - lx*lround((x[i] - x[j])/lx);
		  yij = (y[i] - y[j]) - lx*lround((y[i] - y[j])/lx);
		  zij = (z[i] - z[j]) - lx*lround((z[i] - z[j])/lx);

		  rsq = ((xij*xij) + (yij*yij) + (zij*zij));
		 
		  r = sqrt(rsq);

		//  printf("%f\t%f\n", r, rsq); // ***** Only two values less than 1.5 ... 1.4 and 1 **** some prblm ******** //
		  
		//if(rsq <rc1sq){ 
		    if(rsq <= rc2sq){
		      printf("%lf\n", rsq);	  
		     // printf("%f\t%f\n", r, rsq);
		       UljFT = 1.0/pow(r,12.0);
		       UljST = 1.0/pow(r,6.0);
		       Ulj = 4*(UljFT - UljST + 0.25);   
		    //   printf("%lf\t%lf\n", r, Ulj);
		    }
			
		    Uch = -0.5*30.0*1.5*1.5*log(1-pow((r/1.5), 2.0));
		    //printf("%f\t%f\n", r, Uch);
		  // }
	       }
	    }
	}*/
		
//	printf("%d\n\n", n);
	for(i=1; i<=num;i++){

	   //SAVING TO OLD X,Y,Z
	   xold[i] = x[i];
	   yold[i] = y[i];
	   zold[i] = z[i];
	   
	   //SAVING TO OLD X,Y,Z
	   vxold[i] = vx[i];		      
	   vyold[i] = vy[i];
     	   vzold[i] = vz[i];
    	   //printf("Fe %f\t%f\t%f\t%f\t%f\t%f\n",x[i], y[i], z[i], vxold[i],vyold[i],vzold[i]);
    	  // printf("Fe %f\t%f\t%f\n",x[i], y[i], z[i]);
	   //  printf("Si %f\t%f\t%f\n", x[i] -lx*lround(x[i]/lx), y[i]-lx*lround(y[i]/lx), z[i]-lx*lround(z[i]/lx));
	}

	//-----------APARNA: NOW I KNOW THE FORCE AND POSITION OF MDSTEP 0----------------------------//
	//-----------APARNA: NOW I WILL PROCEED TO CALCULATE THE MDSTEP 1-----------------------------//


	//--------------- Start MD loop --------------------//

	for(mdstep=1; mdstep<=eqmdsteps; mdstep++){

	   //EACH STEP SUM OVER THE PARTICLES AND ITS NEIGBOURS AND
	   //FIND THE FORCE ON EACH PARTICLE AND ENERGY OF THE SYSTEM
	   //AFTER FINDING FORCE FIND NEW POSITION AND VELOCITY
	   for(i=1;i<=num;i++){
	   Fxpar = 0.0;
	   Fypar = 0.0;
     	   Fzpar = 0.0;
	      //SUM OVER NEIGHBOUR IS DONE SO MAKE IT ZERO
	      for(j=1;j<=num;j++){
	         if( i!=j ){
		   //BETTER TO BUILD A NEIGHBOURLIST THAN GOING THROUGH ALL PARTICLES
		   xij = (x[i] - x[j]) - lx*lround((x[i]-x[j])/lx); 
  		   yij = (y[i] - y[j]) - lx*lround((y[i]-y[j])/lx); 
  		   zij = (z[i] - z[j]) - lx*lround((z[i]-z[j])/lx); 
		 
		   rsq = (xij*xij)+(yij*yij)+(zij*zij);
		 
		   // **** Force and energy calculation ****** //
	             r = sqrt(rsq); 

		     if(rsq <= rc2sq){
			//printf("%f\n", r);
		       //Force corresponding to Repulsive LJ potential // 
		       //APARNA: PLEASE CROSS-CHECK AGAIN THE FORCE JUST TO MAKE SURE
		       Forcexlj = 4*((12.0/pow(r,14.0)) - (6.0/pow(r,8.0)))*xij;
		       Forceylj = 4*((12.0/pow(r,14.0)) - (6.0/pow(r,8.0)))*yij;
		       Forcezlj = 4*((12.0/pow(r,14.0)) - (6.0/pow(r,8.0)))*zij;
		       
 		    //	printf("%f\t%f\t%f\t%f\n", r,Forcexlj, Forceylj, Forcezlj); 
		 
		       Fxpar = Fxpar + Forcexlj;
		       Fypar = Fypar + Forceylj;
		       Fzpar = Fzpar + Forcezlj;
		      
		       //F =  (Fxpar* Fxpar) + (Fypar* Fypar) + (Fzpar*Fzpar);
		       //printf("%f\t%f\n", r, sqrt(F));
		     }else{
		       Forcexlj = 0;
		       Forceylj = 0;
		       Forcezlj = 0;
		     }
	
		 }
	    
	      }
	    
	      Fxold[i] = Fxpar; 
   	      Fyold[i] = Fypar;
   	      Fzold[i] = Fzpar;  
	     
	   //  F = (Fxold[i] * Fxold[i]) + (Fyold[i]*Fyold[i]) + (Fzold[i]*Fzold[i]);
	   //   printf("%f\t%f\n", r, sqrt(F));
	   //printf("%f\t%f\t%f\t%f\n", r, Fxold[i], Fyold[i], Fzold[i]);
	   }
	
	// ********** Checked ****** Fine till here ***************** //
	   // *************** Energy and Force calculation done ************************ //
	   // ************* Integrating equation of motion (Velocity verlet) ********** //
	
	   //FIRST FIND THE NEW POSITION
	   for(i=1;i<=num;i++){
              x[i] = xold[i] + vxold[i]*dt + (Fxold[i]/2.0)*dt*dt;                     
	      y[i] = yold[i] + vyold[i]*dt + (Fyold[i]/2.0)*dt*dt;                     
	      z[i] = zold[i] + vzold[i]*dt + (Fzold[i]/2.0)*dt*dt;

	     // printf("Si %f\t%f\t%f\n", x[i], y[i], z[i]);

	      //SET THE NEW TO OLD FOR NEXT LOOP
	      xold[i] = x[i];
	      yold[i] = y[i];
	      zold[i] = z[i];
	   }
	
	   //NOW FIND THE FORCE FROM THE ABOVE NEW POSITION
	   for(i=1;i<=num;i++){
	      //SUM OVER NEIGHBOUR IS DONE SO MAKE IT ZERO
	    Fxpar = 0.0;
	    Fypar = 0.0;
	    Fzpar = 0.0;
	      for(j=1;j<=num;j++){
	         if( i!=j ){
		   //BETTER TO BUILD A NEIGHBOURLIST THAN GOING THROUGH ALL PARTICLES
		   xij = (x[i] - x[j]) - lx*lround((x[i]-x[j])/lx); 
  		   yij = (y[i] - y[j]) - lx*lround((y[i]-y[j])/lx); 
  		   zij = (z[i] - z[j]) - lx*lround((z[i]-z[j])/lx); 
		 
		   rsq = (xij*xij)+(yij*yij)+(zij*zij);
			
		   // **** Force and energy calculation ****** //
	             r = sqrt(rsq);
		     if(rsq <= rc2sq){

		       //Force corresponding to Repulsive LJ potential // 
		       //APARNA: PLEASE CROSS-CHECK AGAIN THE FORCE JUST TO MAKE SURE
		       Forcexlj = 4*((12.0/pow(r,14.0)) - (6.0/pow(r,8.0)))*xij;
		       Forceylj = 4*((12.0/pow(r,14.0)) - (6.0/pow(r,8.0)))*yij;
		       Forcezlj = 4*((12.0/pow(r,14.0)) - (6.0/pow(r,8.0)))*zij;
		       
		       Fxpar = Fxpar + Forcexlj;
		       Fypar = Fypar + Forceylj;
		       Fzpar = Fzpar + Forcezlj;

		       rmin[i] = r;
		      
		       //printf("%f\t%f\t%f\t%f\n", r,Forcexlj, Forceylj, Forcezlj);
		     }else{ //APARNA: IF rsq >= rc2sq, Ulj MIGHT HAVE HOLD SOME PREVIOUS VALUE
			  Forcexlj = 0.0;
			  Forceylj = 0.0;
			  Forcezlj = 0.0;
		     }
		    

		 }
	      }
	   //NEW SET OF FORCES - USED TO CALCULATE THE NEW VELOCITIES
	   Fx[i] = Fxpar;
     	   Fy[i] = Fypar;
	   Fz[i] = Fzpar;
	   //printf("%f\t%f\t%f\n", Fx[i], Fy[i], Fz[i]);
	   }

	   //NOW THE NEW FORCES ARE FOUND... GO ON TO CALCULATE NEW VELOCITIES

	      //printf("%d\n\n", n);
	   for(i=1;i<=num;i++){
	      vx[i] = vxold[i] + ((Fxold[i] + Fx[i])/2.0)*dt;
	      vy[i] = vyold[i] + ((Fyold[i] + Fy[i])/2.0)*dt;
	      vz[i] = vzold[i] + ((Fzold[i] + Fz[i])/2.0)*dt;
	     // printf("%f\t%f\t%f\n", vx[i], vy[i], vz[i]);
	      //SET THE NEW TO OLD FOR NEXT LOOP
	      vxold[i] = vx[i];
	      vyold[i] = vy[i];
	      vzold[i] = vz[i];
	      

	/* begin JAD */
	      //vxcum = 0.0;
	      //vycum = 0.0;
	      //vzcum = 0.0;
	      //
	      //vxcum = vxcum + vx[i];
	      //vycum = vycum + vy[i];
	      //vzcum = vzcum + vz[i];
	/* end JAD */
	   }
	/* begin JAD */
	   //vxcum = vxcum/num;
	   //vycum = vycum/num;
	   //vzcum = vzcum/num;
	   //
	   //for(i=1;i<=num;i++){	   	
	   //   vx[i] = vx[i] - vxcum;
	   //   vy[i] = vy[i] - vycum;
	   //   vz[i] = vz[i] - vzcum;
	   //  //printf("Si %f\t%f\t%f\t%f\t%f\t%f\n", x[i], y[i], z[i], vx[i], vy[i], vz[i]);
	   //  //printf("Si %f\t%f\t%f\n", x[i] -lx*lround(x[i]/lx), y[i]-lx*lround(y[i]/lx), z[i]-lx*lround(z[i]/lx));
	   //}
	/* end JAD */

       // *************** Checked *********** Till here Fine ******************************************* //

	   //THERMO OUTPUT
	   if(mdstep%10 == 0){
	     //printf("Hi\n");
	     //exit(0);
	     //CALCULATE THE POTENTIAL ENERGY AND KINETIC ENERGY
	     ke = 0.0;
	/* begin JAD */
	     pe = 0.0;
	/* end JAD */
	     for(i=1;i<=num;i++){
	/* begin JAD */
	        //pe = 0.0;
	/* end JAD */
		  // printf("Si %f\t%f\t%f\n", x[i], y[i], z[i]);	   
                for(j=1;j<=num;j++){
		   if(i!=j){ 
                     xij = (x[i] - x[j]) - lx*lround((x[i]-x[j])/lx);
		     yij = (y[i] - y[j]) - ly*lround((y[i]-y[j])/ly);
		     zij = (z[i] - z[j]) - lz*lround((z[i]-z[j])/lz);
		    
		     rsq = (xij*xij)+(yij*yij)+(zij*zij);
	              r = sqrt(rsq); 
		      //printf("%f\n", r);
		       //i) Repulsive LJ potential Calculation  //
		       if(rsq <= rc2sq){
                         UljFT = 1.0/pow(r,12.0);
			 UljST = 1.0/pow(r,6.0);
			 Ulj  = 4*(UljFT - UljST + 0.25);
		        //  printf("%f\t%f\n", r, Ulj);
			}
			/* begin JAD */
			else Ulj = 0.0;
			/* end JAD */
		       
		       pe = pe + Ulj ;
		     }
		}
		
	/* begin JAD */
		//ke = ke + 0.5*vx[i]*vx[i];
		ke = ke + 0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	/* end JAD */
		//printf("%f\t%f\t%f\n", vx[i],vy[i],vz[i]);
		//printf("%f\n",ke);
		}

	/* begin JAD */
		pe = pe / 2.0;
	/* end JAD */

	     printf("%d\t%f\t%f\t%f\n", mdstep, pe+ke, pe, ke);
	   }
	}
}

	   //DATA OUTPUT
	 /*  if(mdstep%1000==0){
	   }*/
//	return 0;
//}
