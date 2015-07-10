//Dooley_v3.1
//fix the surrounding molecules around a cavity when calculate the potential
//use the shell model, insertion of water in outer shell are all based on rdf potential, while inner shell are based on actual potential
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
#define  RC_Ueff 10.0 /*unit:A*/
#define  RC_Sel_in 3.0 /*unit:A*/
#define  RC_Sel_out 3.5 /*unit:A*/
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

void copy_array_2(double origin[3][4], double new[3][4])
{
	int x,y;
    for(x=0;x<3;x++)
    	for(y=0;y<4;y++)
    		new[x][y] = origin[x][y];
}

void copy_array_3(double origin[4], double new[4])
{
	int x;
    for(x=0;x<4;x++)
    	new[x] = origin[x];
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

void read_rdf(double g[][4], char addressOO[], char addressOH[], char addressHH[], int nbin)
{
	/*g[distance, OO, OH, HH]*/
	int i;
	char line[100];
	FILE *fp;
	fp = fopen(addressOO, "r");
	for(i=0;i<5;i++)
		fgets(line, 100, fp);
	for(i=0;i<nbin;i++)
	{
		fgets(line, 100, fp);
		sscanf(&line[5],"%lf%lf",&g[i][0],&g[i][1]);
	}
	fclose(fp);

	fp = fopen(addressOH, "r");
	for(i=0;i<5;i++)
		fgets(line, 100, fp);
	for(i=0;i<nbin;i++)
	{
		fgets(line, 100, fp);
		sscanf(&line[10],"%lf",&g[i][2]);
	}
	fclose(fp);

	fp = fopen(addressHH, "r");
	for(i=0;i<5;i++)
		fgets(line, 100, fp);
	for(i=0;i<nbin;i++)
	{
		fgets(line, 100, fp);
		sscanf(&line[10],"%lf",&g[i][3]);
	}
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

void rand_sphere_1(double center[], double coordinate[], double r)  /*generate a vector in a sphere*/
{
	double rsq = 1000.0;
	while(rsq >= r*r)
	{
		coordinate[0] = 2*r*(frand()-0.5);
		coordinate[1] = 2*r*(frand()-0.5);
		coordinate[2] = 2*r*(frand()-0.5);
		rsq = coordinate[0]*coordinate[0]+coordinate[1]*coordinate[1]+coordinate[2]*coordinate[2];
	}
	coordinate[0] += center[0];
	coordinate[1] += center[1];
	coordinate[2] += center[2];
}

void rand_shell(double center[], double coordinate[], double r_in, double r_out)  /*generate a vector in a spherical shell*/
{
	double rsq = 1000.0;
	while(rsq >= r_out*r_out || rsq <= r_in*r_in)
	{
		coordinate[0] = 2*r_out*(frand()-0.5);
		coordinate[1] = 2*r_out*(frand()-0.5);
		coordinate[2] = 2*r_out*(frand()-0.5);
		rsq = coordinate[0]*coordinate[0]+coordinate[1]*coordinate[1]+coordinate[2]*coordinate[2];
	}
	coordinate[0] += center[0];
	coordinate[1] += center[1];
	coordinate[2] += center[2];
}

void rand_sphere_2(double center[], double coordinate[], double r)  /*generate a vector on a sphere surface*/
{
	double rsq = 2.0 ,ran1,ran2,ranh;
	while(rsq >= 1)
	{
		ran1 = 2*(frand()-0.5);
		ran2 = 2*(frand()-0.5);
		rsq = ran1*ran1+ran2*ran2;
	}
	ranh = 2*sqrt(1-rsq);
	coordinate[0] = ran1*ranh*r + center[0];
	coordinate[1] = ran2*ranh*r + center[1];
	coordinate[2] = (1-2*rsq)*r + center[2];
}

void gen_H2(double coord[][4], double L)  /*generate the second H atom in a water molecule in a cavity*/
{
	double a[3], b[3], R11, R12, R13, R21, R22, R23, R31, R32, R33, l = L_OH/10, theta = A_HOH * PI/180, beta, sq;
	distance_pbc(coord[1], coord[0], L, a);
	sq = sqrt(a[1]*a[1]+a[2]*a[2]);
	R11 = a[0]/l;
	R12 = 0;
	R13 = -sq/l;
	R21 = a[1]/l;
	R22 = a[2]/sq;
	R23 = a[0]*a[1]/(l*sq);
	R31 = a[2]/l;
	R32 = -a[1]/sq;
	R33 = a[0]*a[2]/(l*sq);
	beta = (frand()-0.5)*2*PI;
	b[0] = l*cos(theta);
	b[1] = l*sin(theta)*cos(beta);
	b[2] = l*sin(theta)*sin(beta);
	coord[2][0] = b[0]*R11 + b[1]*R12 + b[2]*R13 + coord[0][0];
	coord[2][1] = b[0]*R21 + b[1]*R22 + b[2]*R23 + coord[0][1];
	coord[2][2] = b[0]*R31 + b[1]*R32 + b[2]*R33 + coord[0][2];
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

int serach_number(int n,int array[],int len)       /*used for check a int number is in a array or not,and return the index of the number in the array*/
{
	int i,judge = -1;
	for(i=0;i<len;i++)
	{
		if(n == array[i])
		{
			judge = i;
			break;
		}
	}
	return judge;
}

double Usite_cav(double coordinate[][3][4], int imol, int isite, int N, double L, int index_sel[], int nsel, int isel)           /*used for calculating a particle in cavity*/
{           															/*index_sel takes the index of molecule to be deleted, isel is the order for insertion*/
	int tmol, tsite, s;
	double Ucoulomb=0.0,ULJ=0.0,r,rc;
	for(tmol=0;tmol<N;tmol++)
	{
		s = serach_number(tmol,index_sel,nsel);
		if(tmol==imol)
			continue;
		else if((s!=-1)&&(s>isel))
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
	return ULJ+Ucoulomb;
}

double Urdf(double g[][4], int nbin, double r, int isite, int tsite)  /*input a distance and output corresponding ueff */
{
	int i=0;
		for(i=0;i<nbin;i++)
			if(g[i][0]>r)
				break;
		if(isite + tsite == 0)   /*O-O*/
			return g[i][1];
		else if(isite+tsite>2||(isite==1&&tsite==1))  /*H-H*/
			return g[i][3];
		else    /*O-H*/
			return g[i][2];
}


double Ueff(double coordinate[][3][4], double g[][4], int nbin, int imol, int isite, int N, double L)  /*calculate corresponding Ueff(for one site) based on rdf */
{
	int tmol, tsite;
	double ueff=1.0, r, rc;
	for(tmol=0;tmol<N;tmol++)
	{
		if(tmol==imol)
			continue;
		else
		{
			for(tsite=0;tsite<3;tsite++)
			{
				rc = distance(coordinate[imol][0], coordinate[tmol][0], L);
				if(rc < RC_Ueff/10)
				{
					r = distance(coordinate[imol][isite], coordinate[tmol][tsite], L);
					ueff *=  Urdf(g, nbin, r, isite, tsite);
				}
			}
		}
	}
	return ueff;
}

double Ueff_cav(double coordinate[][3][4], double g[][4], int nbin, int imol, int isite, int N, double L, int index_sel[], int nsel, int isel)
{
	int tmol, tsite, s;
	double ueff=1.0, r, rc;
	for(tmol=0;tmol<N;tmol++)
	{
		s = serach_number(tmol,index_sel,nsel);
		if(tmol==imol)
			continue;
		else if((s!=-1)&&(s>isel))
			continue;
		else
		{
			for(tsite=0;tsite<3;tsite++)
			{
				rc = distance(coordinate[imol][0], coordinate[tmol][0], L);
				if(rc < RC_Ueff/10)
				{
					r = distance(coordinate[imol][isite], coordinate[tmol][tsite], L);
//					printf("g: %.10f\n", Urdf(g, nbin, r, isite, tsite));
					ueff *=  Urdf(g, nbin, r, isite, tsite);
//					printf("product of g: %.10f\n", ueff);
//					if(ueff-0<0.000000000001)
//						for(s=0;s<100;s++)
//							s=1;
				}
			}
		}
	}
	return ueff;
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

int Select(double w[], double sumw)
{
	int n;
	double ws,cumw;
	ws = frand()*sumw;
	cumw = w[0];
	n = 0;
	while(cumw < ws)
	{
		n++;
		cumw += w[n];
	}
	return n;
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
	int N = 4913, nmove = 2000, fout = 20, k = 750,oren = 20, acc_count = 0, natom = 0, num_cps = 0, num_shell = 0, nbin = 1600;
	int n, j,m, jbook, jchoose, nsel, imol=0, isite=0;         /*loop control*/
	int index_sel[N],index_cps[N], index_sel_cps[N];
	double L = 5.27731, T = 300.0, beta, coordinate[N][3][4], temp[N][3][4], center[3],select[k][3][4], selectH[oren][4],box[3],cps[N][3][4],temp_cps[N][3][4], g[nbin][4]; /*cps used to store the molecules in cavity + surrounding*/
	double duration, r, wn,wo,sumw[10],sumwH1,sumwH2,w[10000],wO[10000],wH1[oren],wH2[oren],uH1[oren],uH2[oren],utot[nmove/fout+1],upress[nmove/fout+1];
	double ueff_n,ueff_o,uact_n,uact_o,acc;
	char Title[100], address_out1[100], address_out2[100], address_out3[100], address_out4[100],address_out5[100],address_in[100], address_log[100], addressOO[100], addressOH[100], addressHH[100];
	char directory_in[100] = "/Users/Ziwei/Desktop/INPUT/", open_name[100] ="w100.gro", directory_out[100] = "/Users/Ziwei/Desktop/7_10_1/";
	FILE *fp_log;
	clock_t start, finish;

	start = clock();
	printf("ProName:Dooley_v3.1  Job starts\n");
	n = frand();    /*discard the first random*/
	sprintf(address_log, "%slog.txt", directory_out);
	fp_log = fopen(address_log, "w");
	fprintf(fp_log,"ProName:Dooley_v3.1     JOB STARTS\n");
	fprintf(fp_log,"%d trials will performed on %d molecules at %.2fK\n\n",nmove,N,T);
	fprintf(fp_log,"k = %d, oren = %d, r_cav_in = %.2fA\n, r_cav_out = %.2fA\n",k,oren,RC_Sel_in,RC_Sel_out);
	beta = 1.0/(1.3806488*pow(10,-23)*T);  /*unit for Boltzmann constant is J/K */
	sprintf(address_in, "%s%s", directory_in, open_name);
	fprintf(fp_log,"Open file path:%s\n",address_in);
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
	sprintf(addressOO, "%sstand_rdf_OO.xvg", directory_in);
	sprintf(addressOH, "%sstand_rdf_OH.xvg", directory_in);
	sprintf(addressHH, "%sstand_rdf_HH.xvg", directory_in);
	read_rdf(g, addressOO, addressOH, addressHH, nbin);


	for(n=1;n<(nmove+1);n++)
	{
		nsel = 0;
		num_cps = 0;
		num_shell = 0;
		wn = 1;
		wo = 1;
		center[0] = L*(frand()-0.5);
		center[1] = L*(frand()-0.5);
		center[2] = L*(frand()-0.5);

		for(imol=0;imol<N;imol++)
		{
			r = distance(coordinate[imol][0],center,L);
			if(r < (RC_Sel_out+RC_COULOMB)/10)            /*use RC_Sel+RC_COULOMB to define the cutoff of the surrounding molecules around a cavity*/
			{
				index_cps[num_cps] = imol;
				copy_array_2(coordinate[imol],cps[num_cps]);
				num_cps++;
				if(r < RC_Sel_out/10)
				{
					index_sel[nsel] = imol;
					index_sel_cps[nsel] = num_cps-1;
					nsel++;
					if(r > RC_Sel_in/10)
						num_shell += 1;    /*judge a molecule is inside a shell or outside*/
				}
			}
		}

		copy_array(cps,temp_cps,num_cps);

		for(imol=0;imol<nsel;imol++)
		{
			sumw[imol] = 0;

			for(j=0;j<k;j++)             /*generate k trials for moleclue(oxygen)*/
			{
				if(imol < num_shell)
				{
					w[j] = 1;
					rand_shell(center, cps[index_sel_cps[imol]][0], RC_Sel_in/10, RC_Sel_out/10);
					nojump(L,cps[index_sel_cps[imol]][0]);							/*DECIDE THE POSITION FOR OXYGEN*/
					copy_array_3(cps[index_sel_cps[imol]][0],select[j][0]);
					wO[j] = Ueff_cav(cps, g, nbin, index_sel_cps[imol], 0, num_cps, L, index_sel_cps, nsel,imol);
					w[j] *= wO[j];

					sumwH1 = 0;
					for(m=0;m<oren;m++)
					{
						rand_sphere_2(cps[index_sel_cps[imol]][0], cps[index_sel_cps[imol]][1], L_OH/10);
						nojump(L,cps[index_sel_cps[imol]][1]);							/*DECIDE THE POSITION FOR Hydrogen 1 */
						copy_array_3(cps[index_sel_cps[imol]][1],selectH[m]);
						wH1[m] =  Ueff_cav(cps, g, nbin, index_sel_cps[imol], 1, num_cps, L, index_sel_cps, nsel,imol);
						sumwH1 += wH1[m];
					}
					jbook = Select(wH1,sumwH1);
					copy_array_3(selectH[jbook], select[j][1]);
					copy_array_3(selectH[jbook], cps[index_sel_cps[imol]][1]);
					w[j] *= wH1[jbook];

					sumwH2 = 0;
					for(m=0;m<oren;m++)             /*generate oren trials for hydrogen 2 */
					{
						gen_H2(cps[index_sel_cps[imol]], L);
						nojump(L,cps[index_sel_cps[imol]][2]);							/*DECIDE THE POSITION FOR Hydrogen 2 */
						copy_array_3(cps[index_sel_cps[imol]][2],selectH[m]);
						wH2[m] = Ueff_cav(cps, g, nbin, index_sel_cps[imol], 2, num_cps, L, index_sel_cps, nsel,imol);
						sumwH2 += wH2[m];
					}
					jbook = Select(wH2,sumwH2);
					copy_array_3(selectH[jbook], select[j][2]);
					w[j] *= wH2[jbook];

					sumw[imol] += w[j];
				}
				else
				{
					w[j] = 0;
					rand_sphere_1(center, cps[index_sel_cps[imol]][0], RC_Sel_in/10);
					nojump(L,cps[index_sel_cps[imol]][0]);							/*DECIDE THE POSITION FOR OXYGEN*/
					copy_array_3(cps[index_sel_cps[imol]][0],select[j][0]);
					wO[j] = Usite_cav(cps, index_sel_cps[imol], 0, num_cps, L, index_sel_cps, nsel,imol);
					w[j] += wO[j];

					sumwH1 = 0;
					for(m=0;m<oren;m++)             /*generate oren trials for hydrogen 1 */
					{
						rand_sphere_2(cps[index_sel_cps[imol]][0],  cps[index_sel_cps[imol]][1], L_OH/10);
						nojump(L,cps[index_sel_cps[imol]][1]);							/*DECIDE THE POSITION FOR Hydrogen 1 */
						copy_array_3(cps[index_sel_cps[imol]][1],selectH[m]);
						uH1[m] =  Usite_cav(cps, index_sel_cps[imol], 1, num_cps, L, index_sel_cps, nsel,imol);
						wH1[m] = exp(-beta*uH1[m]*1000/NA);
						sumwH1 += wH1[m];
					}
					jbook = Select(wH1,sumwH1);
					copy_array_3(selectH[jbook], select[j][1]);
					copy_array_3(selectH[jbook], cps[index_sel_cps[imol]][1]);
					w[j] += uH1[jbook];

					sumwH2 = 0;
					for(m=0;m<oren;m++)             /*generate oren trials for hydrogen 2 */
					{
						gen_H2(cps[index_sel_cps[imol]], L);
						nojump(L,cps[index_sel_cps[imol]][2]);							/*DECIDE THE POSITION FOR Hydrogen 2 */
						copy_array_3(cps[index_sel_cps[imol]][2],selectH[m]);
						uH2[m] = Usite_cav(cps, index_sel_cps[imol], 2, num_cps, L, index_sel_cps, nsel,imol);
						wH2[m] = exp(-beta*uH2[m]*1000/NA);
						sumwH2 += wH2[m];
					}
					jbook = Select(wH2,sumwH2);
					copy_array_3(selectH[jbook], select[j][2]);
					w[j] += uH2[jbook];
					w[j] = exp(-beta*w[j]*1000/NA);
					sumw[imol] += w[j];
				}
			}
			jchoose = Select(w,sumw[imol]);
			wn *= sumw[imol];

			sumw[imol] = 0;

			for(j=0;j<k;j++)          /*generate k-1 trials*/
			{
				if(imol < num_shell)
				{
					w[j] = 1;
					if(j != 0)
					{
						rand_shell(center, cps[index_sel_cps[imol]][0], RC_Sel_in/10, RC_Sel_out/10);
						nojump(L,cps[index_sel_cps[imol]][0]);							/*DECIDE THE POSITION FOR OXYGEN*/
						wO[j] = Ueff_cav(cps, g, nbin, index_sel_cps[imol], 0, num_cps, L, index_sel_cps, nsel,imol);
						w[j] *= wO[j];

						sumwH1 = 0;
						for(m=0;m<oren;m++)             /*generate oren trials for hydrogen 1 */
						{
							rand_sphere_2(cps[index_sel_cps[imol]][0], cps[index_sel_cps[imol]][1], L_OH/10);
							nojump(L,cps[index_sel_cps[imol]][1]);							/*DECIDE THE POSITION FOR Hydrogen 1 */
							copy_array_3(cps[index_sel_cps[imol]][1],selectH[m]);
							wH1[m] =  Ueff_cav(cps, g, nbin, index_sel_cps[imol], 1, num_cps, L, index_sel_cps, nsel,imol);
							sumwH1 += wH1[m];
						}
						jbook = Select(wH1,sumwH1);
						copy_array_3(selectH[jbook], cps[index_sel_cps[imol]][1]);

						w[j] *= wH1[jbook];

						sumwH2 = 0;
						for(m=0;m<oren;m++)             /*generate oren trials for hydrogen 2 */
						{
							gen_H2(cps[index_sel_cps[imol]], L);
							nojump(L,cps[index_sel_cps[imol]][2]);							/*DECIDE THE POSITION FOR Hydrogen 2 */
							wH2[m] =  Ueff_cav(cps, g, nbin, index_sel_cps[imol], 2, num_cps, L, index_sel_cps, nsel,imol);
							sumwH2 += wH2[m];
						}
						jbook = Select(wH2,sumwH2);
						w[j] *= wH2[jbook];
					}
					else
					{
						w[j] *=  Ueff_cav(temp_cps, g, nbin, index_sel_cps[imol], 0, num_cps, L, index_sel_cps, nsel,imol);
						w[j] *=  Ueff_cav(temp_cps, g, nbin, index_sel_cps[imol], 1, num_cps, L, index_sel_cps, nsel,imol);
						w[j] *=  Ueff_cav(temp_cps, g, nbin, index_sel_cps[imol], 2, num_cps, L, index_sel_cps, nsel,imol);
					}
				}
				else
				{
					w[j] = 0;
					if(j != 0)
					{
						rand_sphere_1(center, cps[index_sel_cps[imol]][0], RC_Sel_in/10);
						nojump(L,cps[index_sel_cps[imol]][0]);							/*DECIDE THE POSITION FOR OXYGEN*/
						wO[j] = Usite_cav(cps, index_sel_cps[imol], 0, num_cps, L, index_sel_cps, nsel,imol);
						w[j] += wO[j];

						sumwH1 = 0;
						for(m=0;m<oren;m++)             /*generate oren trials for hydrogen 1 */
						{
							rand_sphere_2(cps[index_sel_cps[imol]][0], cps[index_sel_cps[imol]][1], L_OH/10);
							nojump(L,cps[index_sel_cps[imol]][1]);							/*DECIDE THE POSITION FOR Hydrogen 1 */
							copy_array_3(cps[index_sel_cps[imol]][1],selectH[m]);
							uH1[m] = Usite_cav(cps, index_sel_cps[imol], 1, num_cps, L, index_sel_cps, nsel,imol);
							wH1[m] = exp(-beta*uH1[m]*1000/NA);
							sumwH1 += wH1[m];
						}
						jbook = Select(wH1,sumwH1);
						copy_array_3(selectH[jbook], cps[index_sel_cps[imol]][1]);
						w[j] += uH1[jbook];

						sumwH2 = 0;
						for(m=0;m<oren;m++)             /*generate oren trials for hydrogen 2 */
						{
							gen_H2(cps[index_sel_cps[imol]], L);
							nojump(L,cps[index_sel_cps[imol]][2]);							/*DECIDE THE POSITION FOR Hydrogen 2 */
							uH2[m] = Usite_cav(cps, index_sel_cps[imol], 2, num_cps, L, index_sel_cps, nsel,imol);
							wH2[m] = exp(-beta*uH2[m]*1000/NA);
							sumwH2 += wH2[m];
						}
						jbook = Select(wH2,sumwH2);
						w[j] += uH2[jbook];
						w[j] = exp(-beta*w[j]*1000/NA);
					}
					else
					{
						w[j] +=  Usite_cav(temp_cps, index_sel_cps[imol], 0, num_cps, L, index_sel_cps, nsel,imol);
						w[j] +=  Usite_cav(temp_cps, index_sel_cps[imol], 1, num_cps, L, index_sel_cps, nsel,imol);
						w[j] +=  Usite_cav(temp_cps, index_sel_cps[imol], 2, num_cps, L, index_sel_cps, nsel,imol);
						w[j] = exp(-beta*w[j]*1000/NA);
					}
				}
				sumw[imol] += w[j];
			}
			wo *= sumw[imol];

			copy_array_2(select[jchoose], coordinate[index_sel[imol]]);
			copy_array_2(select[jchoose], cps[index_sel_cps[imol]]);
		}

		uact_n = 0;
		uact_o = 0;
		ueff_n = 1;
		ueff_o = 1;
		for(imol=0;imol<num_shell;imol++)
			for(isite=0;isite<3;isite++)
			{
				uact_n += Usite_cav(cps, index_sel_cps[imol], isite, num_cps, L,index_sel_cps, nsel, imol);
				uact_o += Usite_cav(temp_cps, index_sel_cps[imol], isite, num_cps, L,index_sel_cps, nsel, imol);
				ueff_n *= Ueff_cav(cps, g, nbin, index_sel_cps[imol], isite, num_cps, L,index_sel_cps, nsel, imol);
				ueff_o *= Ueff_cav(temp_cps, g, nbin, index_sel_cps[imol], isite, num_cps, L,index_sel_cps, nsel, imol);
			}

		acc = (wn/wo)*(exp(-beta*(uact_n-uact_o)*1000/NA))*(ueff_o/ueff_n);
		/****!!!NOTICE BELOW!!!****/
		if(wn/wo-0<0.0000000001)
			acc = 0;

		printf("wn:%.3f	",wn);
		printf("wo:%.3f	",wo);
		printf("TERM1:%.3f	",(wn/wo));
		printf("uact_n：%.3f	",uact_n);
		printf("uact_o：%.3f	",uact_o);
		printf("TERM2：%.20f	",exp(-beta*(uact_n-uact_o)*1000/NA));
		printf("ueff_o：%.20f	",ueff_o);
		printf("ueff_n：%.20f	",ueff_n);
		printf("TERM3：%.20f	",(ueff_o/ueff_n));
		printf("acc: %.5f\n",acc);


		if(frand() < min(1,acc))
		{
			for(imol=0;imol<nsel;imol++)
				copy_array_2(coordinate[index_sel[imol]], temp[index_sel[imol]]);  /*update temp*/
			acc_count += 1;
			fprintf(fp_log,"ACCEPT\n");
			printf("ACCEPT  ");
		}
		else
		{
			for(imol=0;imol<nsel;imol++)
				copy_array_2(temp[index_sel[imol]], coordinate[index_sel[imol]]); /*set coordinate back to the last one*/
			fprintf(fp_log,"REJECT\n");
			printf("REJECT  ");
		}


		if(n%fout==0)
		{
			sprintf(address_out1, "%swater_out_%d.gro", directory_out,n/fout);
			sprintf(address_out2, "%swater_out_%d.psf", directory_out,n/fout);

//			for(imol=0;imol<nsel;imol++)
//				fprintf(fp_log, "FINALLY!!!  the selected %d molecule:%.3f\n",imol,Usite(coordinate, index_sel[imol], 0, N, L)+Usite(coordinate, index_sel[imol], 1, N, L)+Usite(coordinate, index_sel[imol], 2, N, L));
//				fprintf(fp_log, "FINALLY!!!  the selected %d molecule:%.3f\n",imol,Usite_cav(coordinate, index_sel[imol], 0, N, L, index_sel, nsel,imol)+Usite_cav(coordinate, index_sel[imol], 1, N, L, index_sel, nsel,imol)+Usite_cav(coordinate, index_sel[imol], 2, N, L, index_sel, nsel,imol));

			utot[n/fout] = Utotal(coordinate, N, L);
			upress[n/fout] = Pressure(N, L, beta, coordinate);
			fprintf(fp_log,"The %dth total U is %.3fkJ/mol\n", n/fout, utot[n/fout]);
			fprintf(fp_log,"The %dth system pressure is %.3fatm\n", n/fout, upress[n/fout]);
			output_coord(temp, box, address_out1, address_out2, Title, N*3);
		}
		fprintf(fp_log,"-------------------------------------------------------------------------------------------------------------------------\n");
		printf("The %dth trial is performed, ",n);
		printf("prob %.4f\n",((double)acc_count/(double)n));
	}
	output_utot(address_out3, address_out4, utot, nmove/fout+1, N);
	output_pressure(address_out5, upress, nmove/fout+1);
	output_rdf(coordinate, N, L, directory_out);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	fprintf(fp_log,"\nJob ended, it lasts %.2f seconds, acceptance probability is %.4f\n",duration,((double)acc_count/(double)nmove));
	fclose(fp_log);
}


