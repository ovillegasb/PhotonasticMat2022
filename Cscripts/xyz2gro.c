
/*

Script used to convert STAMP xyz files to GROMACS gro files.

*/

#include <stdio.h>
#include <stdlib.h>

/*

Etapas propuestas:

	1. Leer el archivo xyz.

	2. Extraer informacion del archivo xyz.

	3. Guardar informacion en archivo gro.

*/

int main()
{
	FILE *XYZ;
	FILE *GRO;
	int Natoms;
	char line[256];
	double Lx,Ly,Lz;

	char residue[3] = "RES";

	XYZ = fopen("file.xyz", "r");
	GRO = fopen("file.gro", "w");

	// First line
	// number of atoms
	fgets(line, 256, XYZ);
	Natoms = atoi(line);

	fprintf(GRO, "GRO file, t=0.0\n");
	fprintf(GRO, "%d\n", Natoms);

	printf("number of atoms:\n");
	printf("%d\n", Natoms);

	// Second line
	// Box information
	// Units: angstroms, degree
	fgets(line, 256, XYZ);
	printf("Box information (angstroms):\n");
	sscanf(line, "%lf %lf %lf", &Lx, &Ly, &Lz);
	printf("%f %f %f\n", Lx, Ly, Lz);

	// read the atomic coordinates
	char atsb[5];
	double X, Y, Z;
	for (int i = 1; i <= Natoms; ++i){
		fgets(line, 256, XYZ);
		sscanf(line, "%s %lf %lf %lf", atsb, &X, &Y, &Z);
		fprintf(GRO, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", i, residue, atsb, i, X/10 + Lx/20, Y/10 + Ly/20, Z/10 + Lz/20);
		if (i < 6){
			printf("%d %s %f %f %f\n",i, atsb, X, Y, Z);
		}
	}

	fprintf(GRO, "%10.5f%10.5f%10.5f\n", Lx/10, Ly/10, Lz/10);

	fclose(XYZ);
	fclose(GRO);
	return 0;
}