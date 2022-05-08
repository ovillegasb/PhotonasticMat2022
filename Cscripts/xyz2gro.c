
/*

Script used to convert STAMP xyz files to GROMACS gro files.

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*

Inputs options:
	file.xyz N RES atoms
	argv[1] : file.xyz
	argv[2 3 4] : N1 RES1 atomsRES1
	argv[5 6 7] : N2 RES2 atomsRES2
*/

int main(int argc, char **argv)
{
	FILE *XYZ;
	FILE *GRO;
	int Natoms;
	char line[256];
	double Lx,Ly,Lz;

	// file information
	char name[100];
	char input[100];
	char format[5];
	if (argv[1]==0) {
		printf("No file found");
		exit(0);
	}
	strcpy(input, argv[1]);
	printf("File: %s\n", input);
	char *file = strtok(argv[1], ".");
	strcpy(name, file);
	printf("Name: %s\n", name);

	file = strtok(NULL, ".");
	strcpy(format, file);
	printf("Format: %s\n", format);

	// Residues
	int nres = (argc-2)/3;
	if (nres==0) {
		printf("No residues found");
		exit(0);
	}
	printf("Number of residues: %d\n", nres);

	strcat(name, ".gro");
	printf("Output file: %s\n", name);

	XYZ = fopen(input, "r");
	GRO = fopen(name, "w");

	// First line
	// number of atoms
	fgets(line, 256, XYZ);
	Natoms = atoi(line);

	fprintf(GRO, "GRO file, t=0.0\n");
	fprintf(GRO, "%d\n", Natoms);

	printf("number of atoms in file: %d\n", Natoms);

	// Second line
	// Box information
	// Units: angstroms, degree
	fgets(line, 256, XYZ);
	printf("Box information (angstroms):\n");
	sscanf(line, "%lf %lf %lf", &Lx, &Ly, &Lz);
	printf("%f %f %f\n", Lx, Ly, Lz);

	int natomstot = 0;
	int nmolr;
	char res[3];
	int natomsr;
	for (int n = 0; n < nres; ++n){
		strcpy(res, argv[3*n+3]);
		nmolr = atoi(argv[3*n+2]);
		natomsr = atoi(argv[3*n+4]);
	 	natomstot += natomsr*nmolr;
	}
	printf("Number of atoms introduced: %d\n", natomstot);
	if (natomstot==Natoms){
		printf("The number of atoms entered are equal to the number of atoms in the input file.\n");
	}
	else {
		printf("The number of atoms entered is not equal to the input file.\n");
		exit(0);
	}

	// read the atomic coordinates
	char atsb[5];
	double X, Y, Z;
	int at = 1;
	int atoms_added = 0;
	int ires = 1;
	int ires_added = 0;

	for (int n = 0; n < nres; ++n){
		//printf("%d\n", n);
		strcpy(res, argv[3*n+3]);
		nmolr = atoi(argv[3*n+2]);
		natomsr = atoi(argv[3*n+4]);
		printf("RES: %s\n", res);
	 	printf("NMOL: %d\n", nmolr);
	 	printf("NATOMS: %d\n", natomsr);
	 	printf("NATOMS TOTAL: %d\n", nmolr*natomsr);

	 	while (ires <= nmolr + ires_added){
	 		for (int i = 1; i <= natomsr; ++i){
	 			fgets(line, 256, XYZ);
	 			sscanf(line, "%s %lf %lf %lf", atsb, &X, &Y, &Z);
	 			fprintf(GRO, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", ires, res, atsb, at, X/10 + Lx/20, Y/10 + Ly/20, Z/10 + Lz/20);
	 			at += 1;
	 		}
	 		ires += 1;
	 	}
	 	ires_added += nmolr;
	}

	fprintf(GRO, "%10.5f%10.5f%10.5f\n", Lx/10, Ly/10, Lz/10);

	fclose(XYZ);
	fclose(GRO);
	return 0;
}