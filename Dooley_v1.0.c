//Dooley_v1.0
//conventional MC code, used as reference
//Ziwei Guo  Jul-04-2015


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define  PI 3.1415926
#define  NA 6.022e23
#define  K (8.9875517873681764e9)*(1.6021892e-19)*(1.6021892e-19)*pow(10,9)*NA/1000 /* unit: N·m2/C2 */
#define  EPSILON 0.6364 /*unit:kJ/mol*/
#define  SIGMA 3.15061 /*unit:A*/
#define  L_A 4*EPSILON*pow(SIGMA,12)
#define  L_B 4*EPSILON*pow(SIGMA,6)
#define  RC_COULOMB 12.0  /*unit:A*/
#define  RC_LJ 8.0 /*unit:A*/
#define  ZH +0.417
#define  ZO -0.834
#define  L_OH 0.9572  /*unit:A*/
#define  A_HOH 104.52   /*unit:deg*/
#define  frand() ((double) rand() / (RAND_MAX+1.0))

double min(double a, double b)
	{return ((a)>(b)) ? (b):(a);}

void copy_array(double origin[][3][4], double new[][3][4], int n)
{
	int x,y,z;
	for(x=0;x<n;x++)
    {
    	for(y=0;y<3;y++)
    	{
    		for(z=0;z<4;z++)
    			new[x][y][z] = origin[x][y][z];
    	}
    }
}

void read_coordinate(double coordinate[][3][4], double box[], char address[], char Title[], int natom)
{
	int i;
	char line[100];
	FILE *fp;
	fp = fopen(address, "r");
	fgets(line, 100, fp);
	sscanf(line,"%s", Title);
	fgets(line, 100, fp);
	sscanf(line,"%5d", &natom);
	for (i=0; i<natom; i++)
	{
			fgets(line, 100, fp);
			sscanf(&line[22],"%lf%lf%lf",&coordinate[i/3][i%3][0],&coordinate[i/3][i%3][1],&coordinate[i/3][i%3][2]);
			if(i%3==0)
				coordinate[i/3][i%3][3] = -0.834;
			else
				coordinate[i/3][i%3][3] = +0.417;
	}
	fscanf(fp,"%lf%lf%lf\n",&box[0],&box[1],&box[2]);
	fclose(fp);
}

double distance(double i[], double j[], double L)
{
	double rx = i[0]-j[0], ry = i[1]-j[1], rz = i[2]-j[2];
	if(rx > L/2)
		{rx = rx - L;}
	else if(rx < -L/2)
		{rx = rx + L;}
	if(ry > L/2)
		{ry = ry - L;}
	else if(ry < -L/2)
		{ry = ry + L;}
	if(rz > L/2)
		{rz = rz - L;}
	else if(rz < -L/2)
		{rz = rz + L;}
	return sqrt(rx*rx + ry*ry + rz*rz);
}

int bondjudge_for_psf(double i[], double j[])
{
	int s =1;
	if(fabs(i[0]-j[0])>1 || fabs(i[1]-j[1])>1 || fabs(i[2]-j[2])>1)
		s = 0;
	return s;
}

void distance_pbc(double i[], double j[], double L, double r[])  /*rji  note the order of i,j*/
{
	double rx = i[0]-j[0], ry = i[1]-j[1], rz = i[2]-j[2];
	if(rx > L/2)
		{rx = rx - L;}
	else if(rx < -L/2)
		{rx = rx + L;}
	if(ry > L/2)
		{ry = ry - L;}
	else if(ry < -L/2)
		{ry = ry + L;}
	if(rz > L/2)
		{rz = rz - L;}
	else if(rz < -L/2)
		{rz = rz + L;}
	r[0] = rx;
	r[1] = ry;
	r[2] = rz;
}

void nojump(double L, double r[])  /*ensure this site does not jump out of box !!!Set the origin at the center of the box!!!*/
{
	if(r[0]<-L/2)
		{r[0] += L;}
	else if(r[0]>L/2)
		{r[0] -= L;}
	if(r[1]<-L/2)
		{r[1] += L;}
	else if(r[1]>L/2)
		{r[1] -= L;}
	if(r[2]<-L/2)
		{r[2] += L;}
	else if(r[2]>L/2)
		{r[2] -= L;}
}

double coulomb(double z1, double z2, double r) /* unit: kJ/mol */
{
		return K*z1*z2/r;
}

double LJ(double r) /* unit: kJ/mol */
{
		return L_A*pow(1/(10*r),12)-L_B*pow(1/(10*r),6);
}

void translate(double L, double OS[], double coordinate[][4], double New[][4])
{
	int i;
	double x,y,z,X,Y,Z;
	x = OS[0];
	y = OS[1];
	z = OS[2];
	for(i=0;i<3;i++)
	{
		X = coordinate[i][0]-x;
		Y = coordinate[i][1]-y;
		Z = coordinate[i][2]-z;
		New[i][0] = X;
		New[i][1] = Y;
		New[i][2] = Z;
		nojump(L, New[i]);
	}
}

double cal_angle(double i[], double j[])
{
	double r, ri, rj, theta;
	ri = sqrt(i[0]*i[0]+i[1]*i[1]+i[2]*i[2]);
	rj = sqrt(j[0]*j[0]+j[1]*j[1]+j[2]*j[2]);
	r = (i[0]*j[0]+i[1]*j[1]+i[2]*j[2]);
	theta = acos(r/(ri*rj))*180/PI;
	return theta;
}

void rotate(double OS[], double coordinate[], double theta, double New[])
{
	double R11,R12,R13,R21,R22,R23,R31,R32,R33,x,y,z,r,X,Y,Z;
	theta = theta*PI/180;
	x = OS[0];
	y = OS[1];
	z = OS[2];
	r =sqrt(x*x+y*y+z*z);
	x = x/r;
	y = y/r;
	z = z/r;
	R11 = x*x+(1-x*x)*cos(theta);
	R12 = x*y*(1-cos(theta))-z*sin(theta);
	R13 = x*z*(1-cos(theta))+y*sin(theta);
	R21 = y*x*(1-cos(theta))+z*sin(theta);
	R22 = y*y +(1-y*y)*cos(theta);
	R23 = y*z*(1-cos(theta))-x*sin(theta);
	R31 = z*x*(1-cos(theta))-y*sin(theta);
	R32 = z*y*(1-cos(theta))+x*sin(theta);
	R33 = z*z+(1-z*z)*cos(theta);
	X = coordinate[0]*R11+coordinate[1]*R12+coordinate[2]*R13;
	Y = coordinate[0]*R21+coordinate[1]*R22+coordinate[2]*R23;
	Z = coordinate[0]*R31+coordinate[1]*R32+coordinate[2]*R33;
	New[0] = X;
	New[1] = Y;
	New[2] = Z;
}

/* n,i1,i2 are counters, OS store the vector for translation and reference axis for rotation,OS_T are used to translate
 * molecule back. NewT store the molecule coordinate after translation, NewR store the molecule coordinate after rotation,
 * NewF store the final molecule coordinate after translating back.
 */

void rotate_site(double L, int imol, int isite, double coordinate[][3][4], double theta)
{
	int i1,i2,n;
	double OS[4],OS_T[4],NewT[3][4],NewF[3][4];
	if(isite == 0)  /*if O atom rotate*/
	{
		translate(L, coordinate[imol][1],coordinate[imol],NewT);  /*set the origin at H1 atom */
		OS[0]= NewT[2][0];
		OS[1]= NewT[2][1];
		OS[2]= NewT[2][2];
		rotate(OS, NewT[0], theta, NewT[0]);
		for(n=0;n<3;n++)
			OS_T[n] = -coordinate[imol][1][n];
		translate(L, OS_T,NewT,NewF);
	}
	else   /*if H atom rotate*/
	{
		translate(L, coordinate[imol][0],coordinate[imol],NewT);
		if(isite == 1)  /*if H1 atom rotate*/
		{
			OS[0]= NewT[2][0];
			OS[1]= NewT[2][1];
			OS[2]= NewT[2][2];
			rotate(OS, NewT[1], theta, NewT[1]);
		}
		else  /*if H2 atom rotate*/
		{
			OS[0]= NewT[1][0];
			OS[1]= NewT[1][1];
			OS[2]= NewT[1][2];
			rotate(OS, NewT[2], theta, NewT[2]);
		}
		for(n=0;n<3;n++)
			OS_T[n] = -coordinate[imol][0][n];
		translate(L, OS_T,NewT,NewF);
	}
	for(i1=0;i1<3;i1++)
		for(i2=0;i2<3;i2++)
			coordinate[imol][i1][i2] = NewF[i1][i2];
}

double Usite(double coordinate[][3][4], int imol, int isite, int N, double L)
{
	int tmol, tsite;
	double Ucoulomb=0.0,ULJ=0.0,r,rc;
	for(tmol=0;tmol<N;tmol++)
	{
		if(tmol==imol)
			continue;
		else
		{
			for(tsite=0;tsite<3;tsite++)
			{
				rc = distance(coordinate[imol][0], coordinate[tmol][0], L);
				if(rc < RC_COULOMB/10)
				{
					r = distance(coordinate[imol][isite], coordinate[tmol][tsite], L);
					Ucoulomb += coulomb(coordinate[imol][isite][3], coordinate[tmol][tsite][3], r);
				}
				if(tsite==isite && isite==0 && rc < RC_LJ/10)
					ULJ += LJ(rc);
			}
		}
	}
//	printf("ULJ:%.20f   UC:%.20f   Usite:%.20f\n", ULJ,Ucoulomb,ULJ+Ucoulomb);
	return ULJ+Ucoulomb;
}

double Utotal(double coordinate[][3][4], int N, double L)
{
	int imol, isite;
	double tot = 0.0;
	for(imol=0;imol<N;imol++)
		for(isite=0;isite<3;isite++)
			tot += Usite(coordinate, imol, isite, N, L);
	return tot/2;
}

void output_rdf(double coordinate[][3][4], int N, double L, char directory_address[])
{
	int i,ig, nhis,imol,isite,jmol,jsite;
	double delg = 0.01 ,r,vb,nid_OO,nid_HH,nid_OH;
	double gOO[10000],gHH[10000],gOH[10000];
	char address1[100],address2[100],address3[100];
	FILE *fp1, *fp2, *fp3; /*OO,HH,OH*/

	sprintf(address1,"%srdfOO.xvg",directory_address);
	sprintf(address2,"%srdfHH.xvg",directory_address);
	sprintf(address3,"%srdfOH.xvg",directory_address);
	fp1 = fopen(address1,"w");
	fp2 = fopen(address2,"w");
	fp3 = fopen(address3,"w");
	fprintf(fp1,"@   title\"Radial distribution\"\n");
	fprintf(fp1,"@   xaxis label \"r\"\n");
	fprintf(fp1,"@   yaxis label \"\"\n");
	fprintf(fp1,"@TYPE xy\n");
	fprintf(fp1,"@ subtitle \"O - O\"\n");
	fprintf(fp2,"@   title\"Radial distribution\"\n");
	fprintf(fp2,"@   xaxis label \"r\"\n");
	fprintf(fp2,"@   yaxis label \"\"\n");
	fprintf(fp2,"@TYPE xy\n");
	fprintf(fp2,"@ subtitle \"H - H\"\n");
	fprintf(fp3,"@   title\"Radial distribution\"\n");
	fprintf(fp3,"@   xaxis label \"r\"\n");
	fprintf(fp3,"@   yaxis label \"\"\n");
	fprintf(fp3,"@TYPE xy\n");
	fprintf(fp3,"@ subtitle \"O - H\"\n");

	nhis = (int)((L/2)/delg)+1;
	for(i=0;i<nhis;i++)
	{
		gOO[i] = 0;
		gHH[i] = 0;
		gOH[i] = 0;
	}
	for(imol=0;imol<N-1;imol++)
		for(jmol=imol+1;jmol<N;jmol++)
			for(isite=0;isite<3;isite++)
				for(jsite=0;jsite<3;jsite++)
				{
					r = distance(coordinate[imol][isite],coordinate[jmol][jsite],L);
					if(r<L/2)
					{
						ig = (int)(r/delg);
						if(isite+jsite==0)
							gOO[ig] += 2;
						else if(isite+jsite>2||(isite==1&&jsite==1))
							gHH[ig] += 2;
						else
							gOH[ig] += 2;
					}
				}
	for(i=0;i<nhis;i++)
	{
		r = delg*(i+0.5);
		vb = (pow(i+1,3)-pow(i,3))*pow(delg,3);
		nid_OO = (4./3.)*PI*vb*(N/pow(L,3));
		gOO[i]=gOO[i]/(N*nid_OO);
		fprintf(fp1,"%10.3f%10.4f\n",r,gOO[i]);
		nid_HH = (4./3.)*PI*vb*(2*N/pow(L,3));
		gHH[i]=gHH[i]/(2*N*nid_OO);
		fprintf(fp2,"%10.3f%10.4f\n",r,gHH[i]);
		nid_OH = (4./3.)*PI*vb;
		gOH[i]=gOH[i]/nid_OH*pow(L,3)/(2*N*N)/2;
		fprintf(fp3,"%10.3f%10.4f\n",r,gOH[i]);
	}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
}

void F_vector(int imol, int isite, int jmol, int jsite, double L, double coordinate[][3][4], double f_vector[])
{
	/*f_vector[fx, fy, fz, rx, ry, rz]   fij  */
	double z1,z2,r,rc,mag_c=0.0,mag_lj=0.0,mag,d[3];
	z1 = coordinate[imol][isite][3];
	z2 = coordinate[jmol][jsite][3];

	rc = distance(coordinate[imol][0], coordinate[jmol][0], L);
	r = distance(coordinate[imol][isite], coordinate[jmol][jsite], L);
	if(rc < RC_COULOMB/10)
		mag_c = coulomb(z1, z2, r)/(r*10e-9);
	if(isite==jsite && isite==0 && rc < RC_LJ/10)
		mag_lj = (L_A*pow(1/(10*rc),13)-L_B*pow(1/(10*rc),7))*(10e10);
	distance_pbc(coordinate[jmol][jsite], coordinate[imol][isite], L, d);
	mag = -(mag_lj + mag_c)*1000/NA;
	f_vector[0] = (d[0])*mag/r;      /*unit: N*/
	f_vector[1] = (d[1])*mag/r;
	f_vector[2] = (d[2])*mag/r;
	f_vector[3] = d[0];      /*unit: nm*/
	f_vector[4] = d[1];
	f_vector[5] = d[2];
}

double Pressure(int N, double L, double beta, double coordinate[][3][4])
{
	int imol=0,isite=0,jmol=0,jsite=0;
	double vir = 0.0, fij[6];
	for(imol=0;imol<N-1;imol++)
		for(jmol=imol+1;jmol<N;jmol++)
			for(isite=0;isite<3;isite++)
				for(jsite=0;jsite<3;jsite++)
				{
					F_vector(imol, isite, jmol, jsite, L, coordinate, fij);
					vir += fij[0]*fij[3]+fij[1]*fij[4]+fij[2]*fij[5];     /*fij*rij*/
				}
	vir = vir/3;
	return (N*(10e9)/beta+vir)*(10e9)*(10e9)/(L*L*L)/101325;  /*unit: atm*/
}

void output_coord(double coordinate[][3][4], double box[], char address1[], char address2[], char Title[], int natom)
{
	int i,b=0,t=0,residuenr;
	char residuename[]="WATER";
	char* atomname[]={"OW0","HW1","HW2"};
	double mass[]={16.00,1.00,1.00};
	FILE *fp1, *fp2;  /*fp1 generate a gro file, fp2 generate a psf file*/
	fp1 = fopen(address1, "w");
	fprintf(fp1,"%s\n", Title);
	fprintf(fp1,"%5d\n", natom);
	fp2 = fopen(address2, "w");
	fprintf(fp2,"\n");
	fprintf(fp2,"%8d !NATOM: atoms\n", natom);
	for (i=0; i<natom; i++)
	{
		residuenr = i/3+1;
		fprintf(fp1,"%5d%-5s%5s%5d%8.3lf%8.3lf%8.3lf\n",residuenr,residuename,atomname[i%3],i+1,coordinate[i/3][i%3][0],coordinate[i/3][i%3][1],coordinate[i/3][i%3][2]);
		fprintf(fp2,"%8d U%5d    %-5s  %-5s%-5s%10.6lf%14.4lf           0\n",i+1,residuenr,residuename,atomname[i%3],atomname[i%3],coordinate[i/3][i%3][3],mass[i%3]);
	}
	fprintf(fp1,"%10.5lf%10.5lf%10.5lf\n",box[0],box[1],box[2]);
	fclose(fp1);
	fprintf(fp2,"\n");
	fprintf(fp2,"%8d !NBOND: bonds\n", natom/3*2);
	for (i=0; i<natom/3; i++)
	{
		if(bondjudge_for_psf(coordinate[(t+1)/3][0],coordinate[(t+2)/3][1]))
		{
			fprintf(fp2,"%8d%8d", t+1,t+2);
			b++;
			if(b%4==0)
				fprintf(fp2,"\n");
		}
		if(bondjudge_for_psf(coordinate[(t+1)/3][0],coordinate[(t+2)/3][2]))
		{
			fprintf(fp2,"%8d%8d", t+1,t+3);
			b++;
			if(b%4==0)
			fprintf(fp2,"\n");
		}
		t+=3;
	}
	fprintf(fp2,"\n%d",b); /*b is the number of bonds showed in psf file*/
	fclose(fp2);
}

void output_utot(char address1[], char address2[], double u[],int n, int N)
{
	int i;
	FILE *fp1,*fp2;
	fp1 = fopen(address1, "w");
	fp2 = fopen(address2, "w");
	fprintf(fp1,"@   title\"Total Energy\"\n");
	fprintf(fp1,"@   xaxis label \"Index of output\"\n");
	fprintf(fp1,"@   yaxis label \"Utot(*e4 kJ/mol)\"\n");
	fprintf(fp1,"@TYPE xy\n\n");
	fprintf(fp2,"@   title\"Molecular Energy\"\n");
	fprintf(fp2,"@   xaxis label \"Index of output\"\n");
	fprintf(fp2,"@   yaxis label \"Umol(kJ/mol)\"\n");
	fprintf(fp2,"@TYPE xy\n\n");
	for(i=0;i<n;i++)
	{
		fprintf(fp1,"%10d%10.3f\n",i,u[i]/10000);
		fprintf(fp2,"%10d%10.3f\n",i,u[i]/N);
	}
	fclose(fp1);
	fclose(fp2);
}

void output_pressure(char address[], double p[],int n)
{
	int i;
	FILE *fp;
	fp = fopen(address, "w");
	fprintf(fp,"@   title\"Pressure of System\"\n");
	fprintf(fp,"@   xaxis label \"Index of output\"\n");
	fprintf(fp,"@   yaxis label \"Pressure(atm)\"\n");
	fprintf(fp,"@TYPE xy\n\n");
	for(i=0;i<n;i++)
		fprintf(fp,"%10d%10.3f\n",i,p[i]/10000);
	fclose(fp);
}

/********************************************
**********MAIN PROGRAM STARTS HERE***********
********************************************/

int main(void)
{
	/*N:moleclue number; nmove:number of moves in trajectory; fout:frequency of output calculation
	 * L:box dimension; T:temperature; max_angle:max angle of rotation move(unit:angle! not radial);
	 *  coordinate:[molecule][site][x,y,z,charge] is used to store the  current new structure
	 * temp[number of molecule][site][x,y,z,charge] is used to store the  structure
	 *  for the last move; save_ener[nmove][1] is used to store the total energy*/
	int N = 4913, nmove = 100000000, fout = 1000000, acc_count = 0, natom=0;
	int n, imol=0, isite=0;         /*loop control*/
	double L = 5.27731, T = 300.0, max_angle = 45.0, beta, coordinate[N][3][4], temp[N][3][4], box[3];
	double duration,Usite_o,Usite_n,utot[nmove/fout+1],upress[nmove/fout+1],theta,acc;
	char Title[100], address_out1[100], address_out2[100], address_out3[100], address_out4[100],address_out5[100],address_in[100], address_log[100];
	char directory_in[100] = "/Users/Ziwei/Desktop/", open_name[100] = "ini.gro", directory_out[100] = "/Users/Ziwei/Desktop/TEST/";
	FILE *fp_log;
	clock_t start, finish;

	start = clock();
	printf("ProName:Dooley_v1.0  Job starts\n");
	n = frand();    /*discard the first random*/
	sprintf(address_log, "%slog.txt", directory_out);
	fp_log = fopen(address_log, "w");
	fprintf(fp_log,"ProName:Dooley_v1.0     JOB STARTS\n");
	fprintf(fp_log,"%d trials will performed on %d molecules at %.2fK\n\n",nmove,N,T);
	beta = 1.0/(1.3806488*pow(10,-23)*T);  /*unit for Boltzmann constant is J/K */
	sprintf(address_in, "%s%s", directory_in, open_name);
	read_coordinate(coordinate, box, address_in, Title, natom);
	utot[0] = Utotal(coordinate, N, L);
	upress[0] = Pressure(N, L, beta, coordinate);
	sprintf(address_out3, "%stotal_energy.xvg", directory_out);
	sprintf(address_out4, "%smol_energy.xvg", directory_out);
	sprintf(address_out5, "%spressure.xvg", directory_out);
	for(imol=0;imol<N;imol++)
		for(isite=0;isite<3;isite++)
			nojump(L, coordinate[imol][isite]);
	copy_array(coordinate,temp,N);
	for(n=1;n<(nmove+1);n++)
	{
		imol = (int)(N*frand());
		isite = (int)(3*frand());
		Usite_o = Usite(temp, imol, isite, N, L);
		theta = max_angle*(frand()-0.5)*2;
		rotate_site(L, imol, isite, coordinate, theta);  /*get a new coordinate*/
		nojump(L, coordinate[imol][isite]);
		Usite_n = Usite(coordinate, imol, isite, N, L);
		acc = min(1,exp(-beta*(Usite_n-Usite_o)*1000/NA));
		if(frand()<=acc)
		{
			temp[imol][isite][0] = coordinate[imol][isite][0];  /*update temp*/
			temp[imol][isite][1] = coordinate[imol][isite][1];
			temp[imol][isite][2] = coordinate[imol][isite][2];
			acc_count += 1;
//			fprintf(fp_log,"The %dth trial (rotation angle: %.2f)is performed at molecule %d with this move Accepted\n",n,theta,imol);
		}
		else
		{
			coordinate[imol][isite][0] = temp[imol][isite][0];   /*set coordinate back to the last one*/
			coordinate[imol][isite][1] = temp[imol][isite][1];
			coordinate[imol][isite][2] = temp[imol][isite][2];
//			fprintf(fp_log,"The %dth trial is performed with this move Rejected\n",n);
		}
		if(n%fout==0)
		{
			sprintf(address_out1, "%swater_out_%d.gro", directory_out,n/fout);
			sprintf(address_out2, "%swater_out_%d.psf", directory_out,n/fout);
			utot[n/fout] = Utotal(coordinate, N, L);
			upress[n/fout] = Pressure(N, L, beta, coordinate);
			fprintf(fp_log,"The %dth total U is %.3fkJ/mol\n", n/fout, utot[n/fout]);
			fprintf(fp_log,"The %dth system pressure is %.3fatm\n", n/fout, upress[n/fout]);
			output_coord(temp, box, address_out1, address_out2, Title, N*3);
			if((double)acc_count/(double)n>0.5 && max_angle<90)
				max_angle = max_angle*1.05;
			else
				max_angle = max_angle*0.95;
			fprintf(fp_log,"The %dth max angle is %.3f\n", n/fout, max_angle);
		}
		printf("The %dth trial is performed\n",n);
	}
	output_utot(address_out3, address_out4, utot, nmove/fout+1, N);
	output_pressure(address_out5, upress, nmove/fout+1);
	output_rdf(coordinate, N, L, directory_out);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	fprintf(fp_log,"\nJob ended, it lasts %.2f seconds, acceptance probability is %.4f\n",duration,((double)acc_count/(double)nmove));
	fclose(fp_log);
}

